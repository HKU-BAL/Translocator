//
// Created by Ye Wu on 9/15/2019.
//

#ifndef SNIFFLES_REALIGN_H
#define SNIFFLES_REALIGN_H

#include <vector>
#include <string>
#include <map>
#include "../print/IPrinter.h"
#include "../Alignment.h"
#include "Breakpoint.h"
#include "../bioio.hpp"
#include "../tree/IntervallTree.h"
#include <fstream>

struct hist_str{
    long position;
    int hits;
    std::vector<std::string> names;
};

struct tra_str{
    long position;
    int hits;
    int sameStrand_hits = 0;
    int diffStrand_hits = 0;
    std::map<long, std::vector<long>> map_pos;
    std::vector<std::string> names;

};


struct BpRln{
//    int originalSupport;
    bool denovo;
    bool clipped;
    bool processed;
    bool isSameStrand;
    pair<long, long> coordinate;
    pair<long, long> chr_pos;
    pair<std::string, std::string> chr;
    pair<int, int> chr_idx;
    Breakpoint * bp;
    std::map<std::string,read_str> reads;

    BpRln(bool strand, pair<long, long> coordinates, const RefVector ref, Breakpoint * breakpoint){
        denovo = false;
        this->isSameStrand = strand;
        this->coordinate = coordinates;
        chr_pos.first = IPrinter::calc_pos(coordinate.first, ref, chr_idx.first);
        chr_pos.second = IPrinter::calc_pos(coordinate.second, ref, chr_idx.second);
        chr.first = ref[chr_idx.first].RefName;
        chr.second =ref[chr_idx.second].RefName;

        bp = breakpoint;
    }

    bool operator < (const BpRln& bp) const {
        if (chr.first == bp.chr.first)
            return chr_pos.first < bp.chr_pos.first;
        else return chr_idx.first < bp.chr_idx.first;
    }
};

struct rln_str {
    bool strand;
    bool side;
    pair<string, string> ref;
    pair<long, long> ref_pos;
    pair<long, long> coordinate;
    long diff;
    long read_pos;
    Breakpoint * bp;
};


class ClippedRead {

public:
    long ref_stop;
    long ref_start;
    std::string bases;
    int RefId;

    ClippedRead(){
        this->bases = "";
    }
    void setChrId(int RefId);
    void setBases(std::string bases);
    void set_ref_pos(vector<aln_str> split_events);
    std::vector<rln_str> events_rln;
};

void fix_mismatch_read_pos(vector<differences_str> &event_aln, int i, Alignment * tmp_aln);
void detect_bp_for_realn(Breakpoint  *breakpoint, const RefVector ref, vector<BpRln> &bp_rlns, bool denovo);
void store_tra_pos(vector<tra_str> &positions, read_str read, std::string read_name, bool start);
long get_max_pos(map<long, vector<long>> map_pos, int &max);
long get_max_opp_pos(vector<long> opp_pos, int &max);
bool cal_high_error_side(vector<differences_str> &event_aln, long pos, long distance, Alignment * tmp_aln);
int map_read(Alignment  * tmp_aln, BpRln bp, int distance, const bioio::FastaIndex  index, std::ifstream & fasta);
int map_clipped_read(rln_str event_rln, int distance, const bioio::FastaIndex  index, std::ifstream & fasta, string bases);
void add_realign_read(BpRln bp_realign,  Alignment * tmp_aln);
void add_realign_read(string name, rln_str event_rln);
bool detect_gap(Alignment* tmp_aln, vector<aln_str> split_events, long ref_pos, int distance, RefVector ref);
void detect_clipped_reads_rln(BpRln bp_realign,  Alignment* tmp_aln,
        RefVector ref, std::map<std::string, ClippedRead> &mapClippedRead, ofstream &fasta_out);
void realign_across_read(BpRln bp_realign, vector<differences_str> &event_aln, Alignment* tmp_aln,
                         const bioio::FastaIndex  index, std::ifstream & fasta, RefVector ref, ofstream &fasta_out);
void add_clipped_reads(Alignment *& tmp, std::vector<aln_str> events, short type, RefVector ref, IntervallTree& bst, TNode *&root, long read_id);
long get_ref_lengths(int id, RefVector ref);
bool add_missed_tra(vector<aln_str> split_events, BpRln bpRln, Alignment* tmp_aln, RefVector ref);
void add_high_error_reads(Alignment *& tmp, RefVector ref, IntervallTree& bst, TNode *&root, long read_id) ;
void detect_bps_for_realn(vector<Breakpoint*> breakpoints, vector<Breakpoint*> breakpoints_rln, const RefVector ref, vector<BpRln> &bp_rlns);
void minimap2(string target_path, string query_path, RefVector ref, IntervallTree & bst_rln_denovo,  TNode * &root_rln_denovo);
#endif //SNIFFLES_REALIGN_H

