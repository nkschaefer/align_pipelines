#include <getopt.h>
#include <argp.h>
#include <zlib.h>
#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <htslib/vcf.h>
#include <random>
#include <zlib.h>
#include <time.h>

using std::cout;
using std::endl;
using namespace std;

void read_vcf(string& filename, 
    vector<string>& samples,
    map<int, map<int32_t, int> >& depth_hist,
    map<int, map<int, int> >& qual_hist,
    map<int, pair<int, int> >& het_count,
    map<int, pair<int, int> >& miss_count){

    bcf_hdr_t* bcf_header;
    bcf1_t* bcf_record = bcf_init();
    htsFile* bcf_reader = bcf_open(filename.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR interpreting %s as BCF format.\n", filename.c_str());
        exit(1);
    }
    bcf_header = bcf_hdr_read(bcf_reader);
    int num_samples = bcf_hdr_nsamples(bcf_header);
    for (int i = 0; i < num_samples; ++i){
        samples.push_back(bcf_header->samples[i]);
    }
    
    long int nvar = 0;
    
    int32_t* gts = NULL;
    int n_gts = 0;
    int32_t* gqs = NULL;
    int n_gqs = 0;
    int32_t* dps = NULL;
    int n_dps = 0;

    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        string chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);
        
        // Get 0-based coordinate. Coordinates in BAM will also be 0-based. 
        int pos = bcf_record->pos;
        
        // Load ref/alt alleles and other stuff
        // This puts alleles in bcf_record->d.allele[index]
        // Options for parameter 2:
        /*
        BCF_UN_STR  1       // up to ALT inclusive
        BCF_UN_FLT  2       // up to FILTER
        BCF_UN_INFO 4       // up to INFO
        BCF_UN_SHR  (BCF_UN_STR|BCF_UN_FLT|BCF_UN_INFO) // all shared information
        BCF_UN_FMT  8                           // unpack format and each sample
        BCF_UN_IND  BCF_UN_FMT                  // a synonymo of BCF_UN_FMT
        BCF_UN_ALL (BCF_UN_SHR|BCF_UN_FMT) // everything
        */
        
        bcf_unpack(bcf_record, BCF_UN_STR);
        // want to grab info tags inserted to mark chromosome arm
        //bcf_unpack(bcf_record, BCF_UN_SHR);

        // pass = biallelic, no indels, ref/alt both A/C/G/T
        bool pass = true;
        for (int i = 0; i < bcf_record->n_allele; ++i){
            if (strcmp(bcf_record->d.allele[i], "A") != 0 &&
                strcmp(bcf_record->d.allele[i], "C") != 0 &&
                strcmp(bcf_record->d.allele[i], "G") != 0 && 
                strcmp(bcf_record->d.allele[i], "T") != 0){
                pass = false;
                break;
            }
        }
        if (bcf_record->d.allele[0][0] == bcf_record->d.allele[1][0]){
            pass = false;
        }
        
        if (pass){
            
            // Get all available genotypes.
            //int32_t* gts = NULL;
            //int n_gts = 0;
            int nmiss = 0;
            int num_loaded = bcf_get_genotypes(bcf_header, bcf_record, &gts, &n_gts);
            if (num_loaded <= 0){
                fprintf(stderr, "ERROR loading genotypes at %s %d\n", 
                    chrom.c_str(), bcf_record->pos);
                exit(1);
            }
            
            // Assume ploidy = 2
            int ploidy = 2;
            //int ploidy = n_gts / num_samples; 
            
            // Load depths
            //int32_t* dps = NULL;
            //int n_dps = 0;
            int num_dp_loaded = bcf_get_format_int32(bcf_header, bcf_record, "DP",
                &dps, &n_dps);
            
            // Load genotype qualities
            //int32_t* gqs = NULL;
            //int n_gqs = 0;
            int num_gq_loaded = bcf_get_format_int32(bcf_header, bcf_record, "GQ",
                &gqs, &n_gqs);
            
            for (int i = 0; i < num_samples; ++i){
                if (num_dp_loaded < num_samples || isnan(dps[i]) || dps[i] == bcf_int32_missing){
                    // skip
                }
                else{
                    int32_t dp = dps[i];
                    if (depth_hist.count(i) == 0){
                        map<int32_t, int> m;
                        depth_hist.insert(make_pair(i, m));
                    }        
                    if (depth_hist[i].count(dp) == 0){
                        depth_hist[i].insert(make_pair(dp, 0));
                    }
                    depth_hist[i][dp]++;
                }
                if (num_gq_loaded < num_samples || isnan(gqs[i]) || gqs[i] == bcf_int32_missing){
                    // skip
                }     
                else{
                    int32_t gq = gqs[i];
                    if (qual_hist.count(i) == 0){
                        map<int32_t, int> m;
                        qual_hist.insert(make_pair(i, m));
                    }               
                    if (qual_hist[i].count(gq) == 0){
                        qual_hist[i].insert(make_pair(gq, 0));
                    }
                    qual_hist[i][gq]++;
                }
                int32_t* gtptr = gts + i*ploidy;
                if (miss_count.count(i) == 0){
                    miss_count.insert(make_pair(i, make_pair(0,0)));
                }
                miss_count[i].second++;

                if (bcf_gt_is_missing(gtptr[0])){
                    miss_count[i].first++;
                } 
                else{
                    if (het_count.count(i) == 0){
                        het_count.insert(make_pair(i, make_pair(0,0)));
                    }
                    /*
                    het_count[i].second += 2;
                    if (bcf_gt_allele(gtptr[0]) != 0){
                        het_count[i].first++;
                    }
                    if (bcf_gt_allele(gtptr[1]) != 0){
                        het_count[i].first++;
                    }
                    */
                    het_count[i].second++;
                    if (bcf_gt_allele(gtptr[0]) != bcf_gt_allele(gtptr[1])){
                        het_count[i].first++;
                    }
                }
            }
            //free(gqs);            
            //free(dps);
            //free(gts);
        }
        ++nvar;
        if (nvar % 5000 == 0){
            fprintf(stderr, "Processed %d variants\r", nvar);
        }
    }
    fprintf(stderr, "Processed %ld variants\n", nvar);
    
    free(gqs);
    free(dps);
    free(gts); 

    bcf_hdr_destroy(bcf_header);
    bcf_destroy(bcf_record);
    bcf_close(bcf_reader);
}

float logfac(int x){
    static map<int, float> nfac_cache;
    float nfac = 0;
    if (x > 1){
        if (nfac_cache.count(x) > 0){
            nfac = nfac_cache[x];
        }
        else{
            for (int i = x; i >= 1; --i){
                nfac += log2(i);
            }
            nfac_cache.insert(make_pair(x, nfac));
        }
    }
    return nfac;
}

float fisher_hardy(int nhom0, int nhet, int nhom1){
    int n = nhom0 + nhet + nhom1;
    int np = nhet + 2*nhom1;
    int nq = nhet + 2*nhom0;
    float term1 = logfac(n);
    float term2 = logfac(nhom0) + logfac(nhet) + logfac(nhom1);
    float term3 = nhet + logfac(np) + logfac(nq);
    float term4 = logfac(2*n);
    return (term1-term2) + (term3-term4);
}


void filter_vcf(string& filename, string& outfile, 
    map<int, int32_t>& sample_mindp,
    map<int, int32_t>& sample_maxdp, 
    map<int, set<int> >& pop2idx,
    vector<string>& pops,
    float hardy_thresh,
    int num_samples,
    float frac_missing,
    int min_vq,
    int min_gq,
    string& output_format,
    float sample_prob){
    
    map<int, int> samp2pop;
    if (pop2idx.size() == 0){
        set<int> s;
        pop2idx.insert(make_pair(0, s));
        for (int i = 0; i < num_samples; ++i){
            samp2pop.insert(make_pair(i, 0));
            pop2idx[0].insert(i);
        }
    }
    else{
        for (map<int, set<int> >::iterator pi = pop2idx.begin(); pi != pop2idx.end(); ++pi){
            for (set<int>::iterator pi2 = pi->second.begin(); pi2 != pi->second.end(); ++pi2){
                samp2pop.insert(make_pair(*pi2, pi->first));
            }
        }
    }

    bcf_hdr_t* bcf_header;
    bcf1_t* bcf_record = bcf_init();
    htsFile* bcf_reader = bcf_open(filename.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR interpreting %s as BCF format.\n", filename.c_str());
        exit(1);
    }
    bcf_header = bcf_hdr_read(bcf_reader);
    // mode here can be:
    // w = write VCF
    // wz = write gzVCF
    // wb = write BCF
    string fmtstr = "wz";
    if (output_format == "v"){
        fmtstr = "w";
    }
    else if (output_format == "b"){
        fmtstr = "wb";
    }
    htsFile* hts_out = bcf_open(outfile.c_str(), fmtstr.c_str());
    // Write header to output
    //bcf_hdr_t* out_header = bcf_hdr_init("w");
    bcf_hdr_t* out_header = bcf_hdr_dup(bcf_header);
    int ret = bcf_hdr_write(hts_out, out_header);
    if (ret != 0){
        fprintf(stderr, "ERROR writing bcf header to %s\n", outfile.c_str());
        exit(1);
    }
    long int nvar = 0;
    int32_t* gts = NULL;
    int n_gts = 0;
    int32_t* dps = NULL;
    int n_dps = 0;
    int32_t* gqs = NULL;
    int n_gqs = 0;
    
    int n_hardy_fail = 0;
    int n_depth_fail = 0;
    int n_miss_fail = 0;
    int n_qual_fail = 0;
    int n_gq_fail = 0;
    int n_multi_fail = 0;
    int n_indel_fail = 0;
    
    int nvar_pass = 0;

    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        //string chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);
        bcf_unpack(bcf_record, BCF_UN_ALL);
        ++nvar;   
        bool pass = true;
        for (int i = 0; i < bcf_record->n_allele; ++i){
            if (strcmp(bcf_record->d.allele[i], "A") != 0 &&
                strcmp(bcf_record->d.allele[i], "C") != 0 &&
                strcmp(bcf_record->d.allele[i], "G") != 0 && 
                strcmp(bcf_record->d.allele[i], "T") != 0){
                pass = false;
                n_indel_fail++;
                break;
            }
        }
        if (bcf_record->n_allele != 2){
            pass = false;
            n_multi_fail++;
        }
        if (bcf_record->d.allele[0][0] == bcf_record->d.allele[1][0]){
            pass = false;
        }
        if (min_vq > 0 && bcf_record->qual < min_vq){
            pass = false;
            n_qual_fail++;
        }
        if (pass && sample_prob != -1 && sample_prob != 1.0){
            if ((double)rand() / RAND_MAX > sample_prob){
                pass = false;
            }
        }
        if (pass){
            
            vector<int> pop_hom0;
            vector<int> pop_het;
            vector<int> pop_hom1;
            
            for (map<int, set<int> >::iterator pi = pop2idx.begin(); pi != pop2idx.end(); ++pi){
                pop_hom0.push_back(0);
                pop_het.push_back(0);
                pop_hom1.push_back(0);
            }

            // Get all available genotypes.
            //int32_t* gts = NULL;
            //int n_gts = 0;
            int num_loaded = bcf_get_genotypes(bcf_header, bcf_record, &gts, &n_gts);
            if (num_loaded <= 0){
                fprintf(stderr, "ERROR loading genotypes\n"), 
                exit(1);
            }
            
            // Assume ploidy = 2
            int ploidy = 2;
            
            // Load depths
            //int32_t* dps = NULL;
            //int n_dps = 0;
            int num_dp_loaded = bcf_get_format_int32(bcf_header, bcf_record, "DP",
                &dps, &n_dps);
            
            // Load genotype qualities
            //int32_t* gqs = NULL;
            //int n_gqs = 0;
            int num_gq_loaded = bcf_get_format_int32(bcf_header, bcf_record, "GQ",
                &gqs, &n_gqs);
            
            int n_gt_pass = 0;
            bool hardy_fail = false;
            bool gts_altered = false;
            bool gts_alt_depth = false;
            bool gts_alt_qual = false;

            for (int i = 0; i < num_samples; ++i){
                if (bcf_gt_is_missing(gts[i*ploidy]) || 
                    num_dp_loaded < num_samples || isnan(dps[i]) || dps[i] == bcf_int32_missing){
                    // skip
                }
                else{
                    int32_t dp = dps[i];
                    if (sample_mindp.count(i) > 0 && sample_mindp[i] != -1 && 
                        sample_mindp[i] > dp){
                        // Set to missing.
                        gts[i*ploidy] = bcf_gt_missing;       
                        gts[i*ploidy+1] = bcf_gt_missing;
                        gts_altered = true;
                        gts_alt_depth = true;
                    }
                    else if (sample_maxdp.count(i) > 0 && sample_maxdp[i] != -1 && 
                        sample_maxdp[i] < dp){
                        // Set to missing.
                        gts[i*ploidy] = bcf_gt_missing;   
                        gts[i*ploidy+1] = bcf_gt_missing;
                        gts_altered = true;
                        gts_alt_depth = true;
                    }
                    else{
                        ++n_gt_pass;
                    }        
                }
                if (num_gq_loaded < num_samples || isnan(gqs[i]) || gqs[i] == bcf_int32_missing){
                    // skip
                }     
                else{
                    int32_t gq = gqs[i];
                    if (min_gq > 0 && gq < min_gq){
                        // Set to missing.
                        gts[i*ploidy] = bcf_gt_missing;
                        gts[i*ploidy+1] = bcf_gt_missing;
                        gts_altered = true;
                        gts_alt_qual = true;
                    }
                }
                
                int32_t* gtptr = gts + i*ploidy;
                if (bcf_gt_is_missing(gtptr[0]) || hardy_thresh <= 0){
                    // pass
                } 
                else{
                    int pop = samp2pop[i];
                    int allele1 = bcf_gt_allele(gts[i*ploidy]);
                    int allele2 = bcf_gt_allele(gts[i*ploidy+1]); 
                    
                    if (allele1 + allele2 == 0){
                        pop_hom0[pop]++;
                    }
                    else if (allele1 + allele2 == 1){
                        pop_het[pop]++;
                    }
                    else if (allele1 + allele2 == 2){
                        pop_hom1[pop]++;
                    }
                }
            }
            
            int n_test_hardy = 3; 
            if (n_gt_pass > 0 && samp2pop.size() > 0 && hardy_thresh > 0){
                // test for Hardy-Weinberg equilibrium.
                for (int p = 0; p < pop_hom0.size(); ++p){
                    if (pop_hom0[p] + pop_het[p] + pop_hom1[p] >= n_test_hardy){
                        float fisher_p = fisher_hardy(pop_hom0[p], pop_het[p], pop_hom1[p]);
                        
                        if (pow(2, fisher_p) < hardy_thresh){
                            hardy_fail = true;
                            //fprintf(stderr, "%s\t%d\t%d\t%d\t%f\n", pops[p].c_str(),
                            //    pop_hom0[p], pop_het[p], pop_hom1[p], pow(2, fisher_p));
                            ++n_hardy_fail;
                            break;
                        }
                    }
                }
            } 
            if (((frac_missing == 0 && n_gt_pass == num_samples) || 
                1.0 - (float)n_gt_pass / (float)num_samples < frac_missing) && !hardy_fail){
                if (gts_altered){
                    if (gts_alt_depth){
                        n_depth_fail++;
                    }
                    if (gts_alt_qual){
                        n_gq_fail++;
                    }
                    bcf_update_genotypes(out_header, bcf_record, gts, n_gts); 
                }
                int success = bcf_write(hts_out, out_header, bcf_record);
                ++nvar_pass;
            }
            else if (!hardy_fail){
                // missing fraction failure
                ++n_miss_fail;
            }
            //free(gqs);            
            //free(dps);
            //free(gts);
        }
        if (nvar % 5000 == 0){
            fprintf(stderr, "Processed %d variants\r", nvar);
        }
    }
    fprintf(stderr, "Processed %d variants\n", nvar);
    fprintf(stderr, "%d of %d variants removed - not biallelic\n", n_multi_fail, nvar);
    fprintf(stderr, "%d of %d variants removed - indels\n", n_indel_fail, nvar);
    fprintf(stderr, "%d of %d variants removed for failing variant quality filter\n",
        n_qual_fail, nvar);
    fprintf(stderr, "%d of %d variants removed for failing HWE test\n", n_hardy_fail, nvar);
    fprintf(stderr, "%d of %d variants removed due to fraction of missing genotypes\n", n_miss_fail, nvar);
    fprintf(stderr, "%d of %d kept sites had genotypes removed for failing depth filters\n", n_depth_fail, 
        nvar_pass);
    fprintf(stderr, "%d of %d kept sites had genotypes removed for failing genotype quality filter\n", 
        n_gq_fail, nvar_pass);
    free(gqs);
    free(dps);
    free(gts);

    bcf_hdr_destroy(bcf_header);
    bcf_destroy(bcf_record);
    bcf_close(bcf_reader);
    bcf_close(hts_out);
}

void parse_pops(string& popsfile, 
    vector<string>& samples, 
    vector<string>& pops_uniq,
    map<int, set<int> >& pop2samp){
    
    map<string, int> samp2idx;
    for (int i = 0; i < samples.size(); ++i){
        samp2idx.insert(make_pair(samples[i], i));
    }
    ifstream infile(popsfile.c_str());
    string sample;
    string pop;
    map<string, int> pop2idx;
    
    while (infile >> sample >> pop){
        if (pop2idx.count(pop) == 0){
            pop2idx.insert(make_pair(pop, pops_uniq.size()));
            pops_uniq.push_back(pop);
        }
        if (samp2idx.count(sample) == 0){
            fprintf(stderr, "WARNING: sample %s in pops file not found in VCF.\n", sample.c_str());
        }
        else{
            int pop_idx = pop2idx[pop];
            int samp_idx = samp2idx[sample];
            if (pop2samp.count(pop_idx) == 0){
                set<int> s;
                pop2samp.insert(make_pair(pop_idx, s));
            }
            pop2samp[pop_idx].insert(samp_idx);
        }
    }
}

void help(int code){
    fprintf(stderr, "vcf_filter [OPTIONS]\n");
    fprintf(stderr, "--vcf -v The input VCF file to filter\n");
    fprintf(stderr, "--sample -s To downsample VCF, give a float between 0 and 1 that == probability of keeping \
a record\n");
    fprintf(stderr, "--out -o To filter input VCF, give name of output VCF to create.\n");
    fprintf(stderr, "--output_format -O Format to write output (v = VCF, z = gzVCF, b = BCF). Default = z\n");
    fprintf(stderr, "--depth_low -d Cut genotypes with depths lower than this percentile of per-sample distribution (default 0.005).\n");
    fprintf(stderr, "--depth_high -D Cut genotypes with depths higher than this percentile of per-sample distribution (default 0.995). \
Because of the long right tail of the coverage distribution, this will be taken as this percentile of an exponential distribution \
with mean = (the peak of the coverage distribution)^-1.\n"); 
    fprintf(stderr, "--depth_floor -f Minimum allowable depth across all samples (default = 5)\n");
    fprintf(stderr, "--depth_ceil -c Maximum allowable depth across all samples (default = 1000)\n");
    fprintf(stderr, "--varqual -q Minimum quality of variant to pass filter (default = 50)\n");
    fprintf(stderr, "--gq -Q Minimum quality of genotype to pass filter (default = -1 / NA)\n");
    fprintf(stderr, "--pops -P Provide to calculate within-population Hardy-Weinberg statistics.\n");
    fprintf(stderr, "--hardy -H Cutoff for Hardy-Weinberg p-value to keep a site.\n");
    fprintf(stderr, "--missing -m Allow up to this fraction of input genotypes to be missing at a site. Default = 1.0\n");
}

int main(int argc, char* argv[]){
    // Define arguments 
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"sample", required_argument, 0, 's'},
       {"out", required_argument, 0, 'o'},
       {"depth_low", required_argument, 0, 'd'},
       {"depth_high", required_argument, 0, 'D'},
       {"depth_floor", required_argument, 0, 'f'},
       {"depth_ceil", required_argument, 0, 'c'},
       {"pops", required_argument, 0, 'P'},
       {"hardy", required_argument, 0, 'H'},
       {"missing", required_argument, 0, 'm'},
       {"output_format", required_argument, 0, 'O'},
       {"help", optional_argument, 0, 'h'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    int bufsize = 1048576;
    
    string vcf = "";
    float sample = -1;
    string out;
    //float depth_low = 0.005;
    float depth_low = 0.005;
    float depth_high = 0.995;
    int depth_floor = 5;
    int depth_ceil = 1000;
    int min_vq = 50;
    int min_gq = -1;
    string popsfile = "";
    float hardy_perc = -1;
    float missing = 1.0;
    string output_format = "z";

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "v:s:o:O:p:d:D:f:c:q:Q:H:m:P:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'v':
                vcf = optarg;
                break;
            case 's':
                sample = atof(optarg);
                break;
            case 'o':
                out = optarg;
                break;
            case 'O':
                output_format = optarg;
                break;
            case 'd':
                depth_low = atof(optarg);
                break;
            case 'D':
                depth_high = atof(optarg);
                break;
            case 'f':
                depth_floor = atoi(optarg);
                break;
            case 'c':
                depth_ceil = atoi(optarg);
                break;
            case 'q':
                min_vq = atoi(optarg);
                break;
            case 'Q':
                min_gq = atoi(optarg);
                break;
            case 'P':
                popsfile = optarg;
                break;
            case 'H':
                hardy_perc = atof(optarg);
                break;
            case 'm':
                missing = atof(optarg);
                break;
           case 'h':
                help(0);
                break;
            default:
                help(0);
        }    
    }
    
    if (vcf == ""){
        fprintf(stderr, "ERROR: vcf required\n");
        exit(1);
    }
    if (output_format != "v" && output_format != "z" && output_format != "b"){
        fprintf(stderr, "ERROR: unrecognized output format. Valid options are v, z, or b.\n");
        exit(1);
    }
    if (hardy_perc > 1.0 || depth_low > 1.0 || depth_high > 1.0 || sample > 1.0){
        fprintf(stderr, "--hardy, --depth_low, --depth_high, and --sample all must be a float between 0 and 1. If \
a negative value is passed, it disables the option.\n");
        exit(1);
    } 
    if (out != "" && (depth_low <= 0 && depth_high <= 0 && hardy_perc <= 0)){
        fprintf(stderr, "ERROR: output file name given but no filtering criteria were passed\n");
        exit(1);
    } 
    if (popsfile != "" && hardy_perc <= 0){
        fprintf(stderr, "ERROR: pops file for Hardy-Weinberg testing passed, but no Hardy-Weinberg cutoff given.\n");
        exit(1);
    }
    if (missing < 0.0 || missing > 1.0){
        fprintf(stderr, "ERROR: fraction allowable missing genotypes must be between 0 and 1.\n"); 
    }
    if (depth_low >= depth_high){
        fprintf(stderr, "ERROR: lower depth cutoff is higher than upper depth cutoff\n");
        exit(1);   
    }

    // Initialize random number seed    
    srand(time(NULL));
    
    string outbase;
    if (out.rfind(".vcf") != string::npos){
        outbase = out.substr(0, out.rfind(".vcf"));   
    }
    else{
        outbase = out;
        if (output_format == "v"){
            out += ".vcf";
        }
        else if (output_format == "z"){
            out += ".vcf.gz";
        }
        else if (output_format == "b"){
            out += ".bcf";
        }
    }

    srand(time(NULL)); 
    
    vector<string> samples;
    map<int, map<int32_t, int> > depth_hist;
    map<int, map<int, int> > qual_hist;
    map<int, pair<int, int> > het_count;
    map<int, pair<int, int> > miss_count;
    
    read_vcf(vcf, samples, depth_hist, qual_hist, het_count, miss_count);
    
    map<int, set<int> > pop2idx;
    vector<string> popnames;
    bool pops = false;
    if (popsfile != ""){
        pops = true;
        parse_pops(popsfile, samples, popnames, pop2idx);
    }

    // Print results
    vector<float> sample_depth_sum;
    vector<float> sample_gq_sum;
    
    // prevent overflow
    float sumscale = 1e-6;
    
    string statsname = outbase + ".hist";
    FILE* statsf = fopen(statsname.c_str(), "w");
    
    // Write parameters to stats file 
    fprintf(statsf, "# vcf %s\n", vcf.c_str());
    fprintf(statsf, "# sample %f\n", sample);
    fprintf(statsf, "# depth %f %f\n", depth_low, depth_high);
    fprintf(statsf, "# depth_floor_ceil %d %d\n", depth_floor, depth_ceil);
    fprintf(statsf, "# min_gq %d min_vq %d\n", min_gq, min_vq);
    fprintf(statsf, "# hardy_perc %f pops %s\n", hardy_perc, popsfile.c_str());
    fprintf(statsf, "# frac_missing %f\n", missing);

    for (int i = 0; i < samples.size(); ++i){
        float ds = 0;
        for (map<int32_t, int>::iterator bin = depth_hist[i].begin(); bin != depth_hist[i].end(); ++bin){
            fprintf(statsf, "DP\t%s\t%d\t%d\n", samples[i].c_str(),
                bin->first, bin->second);
            ds += (float)(bin->first * bin->second)/sumscale;
        }
        sample_depth_sum.push_back(ds);
        float gs = 0;
        for (map<int, int>::iterator bin = qual_hist[i].begin(); bin != qual_hist[i].end(); ++bin){
            fprintf(statsf, "GQ\t%s\t%d\t%d\n", samples[i].c_str(), bin->first, bin->second);
            gs += (float)(bin->first * bin->second)/sumscale;
        }
        sample_gq_sum.push_back(gs);
        fprintf(statsf, "HET\t%s\t%d\t%d\n", samples[i].c_str(), het_count[i].first, het_count[i].second);
        fprintf(statsf, "MISS\t%s\t%d\t%d\n", samples[i].c_str(), miss_count[i].first, miss_count[i].second);
    }   
    
    map<int, int32_t> sample_mindp;
    map<int, int32_t> sample_maxdp;

    if (depth_low > 0 || depth_high > 0 || depth_floor > 0 || depth_ceil > 0){
        fprintf(stderr, "Determining sample-specific depth cutoffs...\n");
        
        for (int i = 0; i < samples.size(); ++i){
            float tot = sample_depth_sum[i];
            float min = depth_low * tot;
            float max = depth_high * tot;
            float med = 0.5 * tot;
            int mindp = -1;
            int maxdp = -1;
            int med_dp = -1;
            float cumulative = 0;
            int32_t mode = -1;
            int modecount = -1;
            for (map<int32_t, int>::iterator bin = depth_hist[i].begin(); bin != depth_hist[i].end(); ++bin){
                float next = cumulative + (float)(bin->first * bin->second)/sumscale;
                if (depth_low > 0 && mindp == -1 && next >= min){
                    mindp = bin->first;
                }
                if (med_dp == -1 && next >= med){
                    med_dp = bin->first;
                }
                if (depth_high > 0 && maxdp == -1 && next >= max){
                    maxdp = bin->first-1;
                }
                if (bin->first != 0 && (mode == -1 || bin->second > modecount)){
                    modecount = bin->second;
                    mode = bin->first;         
                }
                cumulative = next;
            }
            
            // New strategy: use empirical lower threshold, but avoid this for
            // upper because of extremely long right tail that is difficult to 
            // model and very data set-specific.
            // Instead, take mode from histogram as 1/mean of exponential distribution
            // and calculate percentile of this distribution as upper cutoff.
            float expmean = 1.0 / (float)mode;
            mindp = (int)round(-log(1.0 - depth_low)/expmean);
            maxdp = (int)round(-log(1.0 - depth_high)/expmean);

            if (depth_floor > 0 && (mindp == -1 || depth_floor > mindp)){
                mindp = depth_floor;
            }
            if (depth_ceil > 0 && (maxdp == -1 || depth_ceil < maxdp)){
                maxdp = depth_ceil;
            }
            if (maxdp <= mindp){
                maxdp = mindp + 1;
            }
            sample_mindp.insert(make_pair(i, mindp));
            sample_maxdp.insert(make_pair(i, maxdp));
            fprintf(statsf, "DP_CUTOFF\t%s\t%d\t%d\n", samples[i].c_str(),
                mindp, maxdp);

            // Obtain median GQ too
            float med2 = 0.5 * sample_gq_sum[i];
            int med_gq = -1;
            float cumulative_gq = 0;
            int32_t mode_gq = -1;
            int mode_gq_count = -1;
            for (map<int32_t, int>::iterator bin = qual_hist[i].begin(); bin != qual_hist[i].end(); ++bin){
                float next = cumulative_gq + (float)(bin->first * bin->second)/sumscale;
                if (med_gq == -1 && next >= med2){
                    med_gq = bin->first;
                }
                if (mode_gq == -1 || bin->second > mode_gq_count){
                    mode_gq_count = bin->second;
                    mode_gq = bin->first;
                }
                cumulative_gq = next;
            }
            if (med_gq != -1 || med_dp != -1){
                fprintf(statsf, "MED_DP_GQ\t%s\t%d\t%d\n", samples[i].c_str(),
                    med_dp, med_gq);  
            }
            if (mode != -1 || mode_gq != -1){
                fprintf(statsf, "MODE_DP_GQ\t%s\t%d\t%d\n", samples[i].c_str(),
                    mode, mode_gq);
            }
        }       
    }
    fclose(statsf);

    depth_hist.clear();
    qual_hist.clear();
    het_count.clear();
    miss_count.clear();
    sample_depth_sum.clear();
    sample_gq_sum.clear();

    if (depth_low > 0 || depth_high > 0 || hardy_perc > 0){
        // Filter input VCF.
        if (out == ""){
            out = "-";
        }
        fprintf(stderr, "Filtering variants...\n");
        filter_vcf(vcf, out, sample_mindp, sample_maxdp, pop2idx, 
            popnames, hardy_perc, samples.size(), missing,
            min_vq, min_gq, output_format, sample);
    } 
}
