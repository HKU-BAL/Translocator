//
// Created by Ye Wu on 9/15/2019.
//

#include "Realign.h"
#include "limits.h"
#include "seqan/align.h"
#include "seqan/bam_io.h"
#include "seqan/modifier.h"
#include <set>
#include <iostream>

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "../../lib/minimap2/minimap.h"
#include "../../lib/minimap2/kseq.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;

void ClippedRead::setChrId(int RefId){
    this->RefId = RefId;
}
void ClippedRead::setBases(string bases){
    this->bases = bases;
}
void ClippedRead::set_ref_pos(vector<aln_str> split_events){
    long max_pos=0, min_pos=LONG_MAX;
    for (aln_str split_event: split_events){
        if (split_event.RefID != this->RefId)
            continue;
        if (split_event.pos < min_pos)
            min_pos = split_event.pos;
        if (split_event.pos + split_event.length > max_pos)
            max_pos = split_event.pos + split_event.length;
    }
    this->ref_start = min_pos;
    this->ref_stop = max_pos;
}

long get_ref_lengths(int id, RefVector ref) {
    long length = 0;

    for (size_t i = 0; i < (size_t) id && i < ref.size(); i++) {
        length += (long) ref[i].RefLength + (long) Parameter::Instance()->max_dist;
    }
    return length;
}

void fix_mismatch_read_pos(vector<differences_str> &event_aln, int i, Alignment * tmp_aln){
    int read_pos = 0;
    if (i == 0) {
        if (tmp_aln->getAlignment()->CigarData[0].Type == 'S')
            read_pos += tmp_aln->getAlignment()->CigarData[0].Length;
        read_pos += event_aln[i].position - tmp_aln->getPosition();
        event_aln[i].readposition = read_pos;
    } else
        event_aln[i].readposition = event_aln[i-1].readposition + event_aln[i].position - event_aln[i-1].position -
                                    event_aln[i-1].type;
}

void detect_bp_for_realn(Breakpoint  *breakpoint, const RefVector ref, vector<BpRln> &bp_rlns, bool denovo) {
    auto point = breakpoint->get_coordinates();
    vector<tra_str> pos_start;

    vector<tra_str> pos_stop;


    for (auto i = point.support.begin(); i != point.support.end(); ++i) {
        store_tra_pos(pos_start, (*i).second, (*i).first, true);
        store_tra_pos(pos_stop, (*i).second, (*i).first, false);
    }

    for (size_t i = 0; i < pos_start.size(); i++) {

        if (pos_start[i].hits >= Parameter::Instance()->min_support / 3) {
//        if (pos_start[i].hits >= 2) {

            bool isSameStrand = pos_start[i].sameStrand_hits >= pos_start[i].diffStrand_hits;
//            bool isSameStrand = true;
            pair<long, long> coordinate;
            int max_start = 0, max_stop = 0;
            coordinate.first = get_max_pos(pos_start[i].map_pos, max_start);
            coordinate.second = get_max_opp_pos(pos_start[i].map_pos[coordinate.first], max_stop);

            BpRln start_bp(isSameStrand, coordinate, ref, breakpoint);
            std::swap(coordinate.first, coordinate.second);
            BpRln stop_bp(isSameStrand, coordinate, ref, breakpoint);
            if (denovo) {
                start_bp.denovo = true;
                bp_rlns.push_back(start_bp);
            } else {
                bp_rlns.push_back(start_bp);
                bp_rlns.push_back(stop_bp);
            }
        }
    }
}

void split(const std::string &s, char delim, vector<string> &result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) result.push_back(item);
    }
}

void store_tra_pos(vector<tra_str> &positions, read_str read, std::string read_name, bool start) {
    long pos, opp_pos;
    if (start) {
        pos = read.coordinates.first;
        opp_pos = read.coordinates.second;
    } else {
        opp_pos = read.coordinates.second;
        pos = read.coordinates.first;
    }


    for (size_t i = 0; i < positions.size(); i++) {
        if (abs(positions[i].position - pos) < Parameter::Instance()->min_length) {
            positions[i].hits++;
            positions[i].names.push_back(read_name);
            positions[i].map_pos[pos].push_back(opp_pos);
            if (read.read_strand.first == read.read_strand.second)
                positions[i].sameStrand_hits++;
            else positions[i].diffStrand_hits++;
            return;
        }
    }
    tra_str tmp;
    tmp.position = pos;
    tmp.hits = 1;
    tmp.names.push_back(read_name);
    tmp.map_pos[pos].push_back(opp_pos);
    if (read.read_strand.first == read.read_strand.second)
        tmp.sameStrand_hits = 1;
    else tmp.diffStrand_hits = 1;
    positions.push_back(tmp);
}

long get_max_pos(map<long, vector<long>> map_pos, int &max){
    long coordinate = 0;
    for (auto i = map_pos.begin(); i != map_pos.end(); i++){
        if ((*i).second.size() > max){
            max = (*i).second.size();
            coordinate = (*i).first;
        }
    }
    return coordinate;
}

long get_max_opp_pos(vector<long> opp_pos, int &max){
    map<long, int> map_opp_pos;
    for (auto j = opp_pos.begin(); j != opp_pos.end(); j++)
        if (map_opp_pos.find(*j) == map_opp_pos.end())
            map_opp_pos[*j] = 1;
        else map_opp_pos[*j]++;
    long coordinate = 0;

    for (auto i = map_opp_pos.begin(); i != map_opp_pos.end(); i++){
        if ((*i).second > max){
            max = (*i).second;
            coordinate = (*i).first;
        }
    }
    return coordinate;
}

bool cal_high_error_side(vector<differences_str> &event_aln, long pos, long distance, Alignment * tmp_aln) {
    if (event_aln.empty())
        return false;
    for (size_t j = 0; j < event_aln.size(); j++) {
        if (event_aln[j].type == 0)
            fix_mismatch_read_pos(event_aln, j, tmp_aln);
    }
    int left_hits = 0;
    int right_hits = 0;
    differences_str evt_left = event_aln[0];
    differences_str evt_right = event_aln[event_aln.size()-1];
    for (auto i : event_aln){
        if (i.position < pos - distance) continue;
        if (i.position > pos + distance) break;
        if (i.position >= pos - distance && i.position <= pos){
            left_hits += max((int) abs(i.type), 1);
            evt_left = i;
        }
        else if (i.position >= pos && i.position <= pos + distance) {
            right_hits += max((int) abs(i.type), 1);
            if (evt_right.position > i.position) evt_right = i;
        }

    }

    if (abs(evt_left.position - pos) < abs(evt_right.position - pos)) {
        tmp_aln->bp_read_pos = evt_left.readposition;
        tmp_aln->diff = evt_left.position - pos;
    } else if (abs(evt_right.position - pos) < abs(evt_left.position - pos)) {
        tmp_aln->bp_read_pos = evt_right.readposition;
        tmp_aln->diff = evt_right.position - pos;
    }

    if (right_hits - left_hits >= distance / 5) {
        tmp_aln->high_error_region_score = -1 * right_hits;
        if (abs(evt_left.position - pos) == abs(evt_right.position - pos)) {
            tmp_aln->bp_read_pos = evt_left.readposition;
            tmp_aln->diff = evt_left.position - pos;
        }
        tmp_aln->high_error_side = true;
        return true;
    }

    if (left_hits - right_hits >= distance / 5) {
        tmp_aln->high_error_region_score = -1 * left_hits;
        if (abs(evt_left.position - pos) == abs(evt_right.position - pos)) {
            tmp_aln->bp_read_pos = evt_right.readposition;
            tmp_aln->diff = evt_right.position - pos;
        }
        tmp_aln->high_error_side = false;
        return true;
    }

    return false;
}

int map_read(Alignment  * tmp_aln, BpRln bp, int distance,
             const bioio::FastaIndex  index, std::ifstream & fasta){
    using namespace seqan;
    typedef String<Dna> TSequence;
    typedef StringSet<TSequence> TStringSet;
    typedef Align<TSequence, ArrayGaps> TAlign;// container for strings
    typedef StringSet<TSequence, Dependent<>> TDepStringSet;
    typedef seqan::Alignment<TDepStringSet> TAlignStringSet; // dependent string set
    typedef Graph<TAlignStringSet> TAlignGraph;       // alignment graph// sequence type


    if (bp.isSameStrand) bp.chr_pos.second += tmp_aln->diff;
    else bp.chr_pos.second  -= tmp_aln->diff;

    if (!tmp_aln->high_error_side) tmp_aln->bp_read_pos -= distance; //lefthand side
    if (bp.isSameStrand != tmp_aln->high_error_side) bp.chr_pos.second -= distance;

    string ref_str = bioio::read_fasta_contig(fasta, index.at(bp.chr.second), bp.chr_pos.second, distance);
    transform(ref_str.begin(), ref_str.end(), ref_str.begin(), ::toupper);
    TSequence reference = ref_str;
    if (tmp_aln->bp_read_pos + distance > tmp_aln->getQueryBases().size())
        return -50;
    if (tmp_aln->bp_read_pos < 0)
        return -50;
//    std::cout << bp.chr.second << " " << bp.chr_pos.second << endl;
//    std::cout << tmp_aln->bp_read_pos << " " << distance << endl;
    string seq_str = tmp_aln->getQueryBases().substr(tmp_aln->bp_read_pos, distance);
    transform(seq_str.begin(), seq_str.end(), seq_str.begin(), ::toupper);
    TSequence sequence = seq_str;
    TStringSet sequences;
    appendValue(sequences, reference);
    appendValue(sequences, sequence);

    TAlignGraph alignG(sequences);

    int score = globalAlignment(alignG, Score<int, Simple>(0, -1, -1), AlignConfig<false, false, true, true>(), LinearGaps());

    return score;

}


int map_clipped_read(rln_str event_rln, int distance, const bioio::FastaIndex  index, std::ifstream & fasta, string bases) {
    using namespace seqan;
    typedef String<char> TSequence;
    typedef StringSet<TSequence> TStringSet;
    typedef Align<TSequence, ArrayGaps> TAlign;// container for strings
    typedef StringSet<TSequence, Dependent<>> TDepStringSet;
    typedef seqan::Alignment<TDepStringSet> TAlignStringSet; // dependent string set
    typedef Graph<TAlignStringSet> TAlignGraph;       // alignment graph// sequence type

    if (event_rln.strand) event_rln.ref_pos.second += event_rln.diff;
    else event_rln.ref_pos.second -= event_rln.diff;

    if (!event_rln.side) event_rln.read_pos -= distance; //lefthand side
    if (event_rln.strand != event_rln.side) event_rln.ref_pos.second -= distance;

    string ref_str = bioio::read_fasta_contig(fasta, index.at(event_rln.ref.second), event_rln.ref_pos.second, distance);
    transform(ref_str.begin(), ref_str.end(), ref_str.begin(), ::toupper);
    TSequence reference = ref_str;
    if (event_rln.read_pos + distance > bases.size())
        return -50;

    if (event_rln.read_pos < 0)
        return -50;
//    std::cout << bp.chr.second << " " << bp.chr_pos.second << endl;
//    std::cout << tmp_aln->bp_read_pos << " " << distance << endl;
    string seq_str = bases.substr(event_rln.read_pos, distance);
    transform(seq_str.begin(), seq_str.end(), seq_str.begin(), ::toupper);
    TSequence sequence = seq_str;
    TStringSet sequences;
//    cout << "seq: " << sequence << endl;
    if (!event_rln.strand)
        reverseComplement(sequence);
    appendValue(sequences, reference);
    appendValue(sequences, sequence);

    TAlignGraph alignG(sequences);

//    cout << "ref: " << reference << endl;
//    cout << event_rln.ref.second << " " << event_rln.ref_pos.second << endl;
//
//    cout << "seq: " << sequence << endl;


    int score = globalAlignment(alignG, Score<int, Simple>(0, -1, -1), AlignConfig<false, false, true, true>(), LinearGaps());

//    std::cout << score << endl;
//    std::cout << alignG << endl;


    return score;


}

void add_realign_read(BpRln bp_realign,  Alignment * tmp_aln){
    read_str tmp;
    if (bp_realign.coordinate.first < bp_realign.coordinate.second) {
        tmp.coordinates.first = bp_realign.coordinate.first + tmp_aln->diff;
        if (bp_realign.isSameStrand)
            tmp.coordinates.second = bp_realign.coordinate.second + tmp_aln->diff;
        else tmp.coordinates.second = bp_realign.coordinate.second - tmp_aln->diff;
    } else {
        tmp.coordinates.second = bp_realign.coordinate.first + tmp_aln->diff;
        if (bp_realign.isSameStrand)
            tmp.coordinates.first = bp_realign.coordinate.second + tmp_aln->diff;
        else tmp.coordinates.first = bp_realign.coordinate.second - tmp_aln->diff;
    }
    tmp.SV = TRA;
    tmp.type = 1;
//                        std::cout << bp_realign.bp->get_coordinates().support.size() << endl;
    auto map = bp_realign.bp->get_coordinates().support;
    if (map.find(tmp_aln->getName()+"_rln_0") == map.end())
        bp_realign.bp->add_read(tmp, tmp_aln->getName()+"_rln_0");
    else
        bp_realign.bp->add_read(tmp, tmp_aln->getName() + "_rln_1");
}

void add_realign_read(string name, rln_str event_rln){
    read_str tmp;
    if (event_rln.coordinate.first < event_rln.coordinate.second) {
        tmp.coordinates.first = event_rln.coordinate.first + event_rln.diff;
        tmp.coordinates.second = event_rln.coordinate.second + event_rln.diff;
    } else {
        tmp.coordinates.second = event_rln.coordinate.first + event_rln.diff;
        if (event_rln.strand)
            tmp.coordinates.first = event_rln.coordinate.second + event_rln.diff;
        else tmp.coordinates.first = event_rln.coordinate.second - event_rln.diff;
    }
    tmp.SV = TRA;
    tmp.type = 1;

//                        std::cout << bp_realign.bp->get_coordinates().support.size() << endl;
    auto map = event_rln.bp->get_coordinates().support;
    if (map.find(name+"_rln_0") == map.end())
        event_rln.bp->add_read(tmp, name+"_rln_0");
    else
        event_rln.bp->add_read(tmp, name + "_rln_1");
}


bool detect_gap(Alignment* tmp_aln, std::vector <aln_str> split_events, long ref_pos, int distance, RefVector ref) {

//    cout << ref_pos << endl;
//    cout << "step2" << endl;
    if (split_events.size() == 0) return false;
    long read_length;
    aln_str last_event = split_events[split_events.size()-1];
//    cout << "step2.1" << endl;
    if (tmp_aln->getAlignment()->IsPrimaryAlignment() && !(tmp_aln->getAlignment()->AlignmentFlag & 0x800) )
        read_length = tmp_aln->getAlignment()->Length;
    else if (last_event.strand) {
        read_length = last_event.read_pos_stop;
        if (last_event.cigar.size() > 0 && last_event.cigar.back().Type == 'S')
            read_length += last_event.cigar.back().Length;
    } else {
        read_length = last_event.read_pos_stop;
        if (last_event.cigar.size() > 0 && last_event.cigar.front().Type == 'S')
            read_length += last_event.cigar.front().Length;
    }
//    cout << "step2.2" << endl;
    if (split_events[0].isMain && split_events[0].read_pos_start >= distance) {
        if (split_events[0].strand && abs(split_events[0].pos - ref_pos) < distance) {
            tmp_aln->diff = split_events[0].pos - ref_pos;
            tmp_aln->high_error_side = false;
            tmp_aln->bp_read_pos = split_events[0].read_pos_start;
            return true;
        } else if (!split_events[0].strand && abs(split_events[0].pos + split_events[0].length - ref_pos) < distance) {
            tmp_aln->diff = split_events[0].pos + split_events[0].length - ref_pos;
            tmp_aln->high_error_side = true;
            tmp_aln->bp_read_pos = read_length - 1 - split_events[0].read_pos_start;
            return true;
        }
    }

    if (last_event.isMain && tmp_aln->getAlignment()->Length - 1 - last_event.read_pos_stop >= distance) {
        if (last_event.strand && abs(last_event.pos + last_event.length - ref_pos) < distance) {
            tmp_aln->diff = last_event.pos + last_event.length - ref_pos;
            tmp_aln->high_error_side = true;
            tmp_aln->bp_read_pos = last_event.read_pos_stop;
            return true;
        } else if (!last_event.strand && abs(last_event.pos - ref_pos) < distance) {
            tmp_aln->diff = last_event.pos - ref_pos;
            tmp_aln->high_error_side = false;
            tmp_aln->bp_read_pos = read_length - 1 - last_event.read_pos_stop;
            return true;
        }
    }
//    cout << "step2.2" << endl;
    for (size_t i = 0; i < split_events.size(); i++) {
        if (split_events[i].isMain) {
            if (abs(split_events[i].pos - ref_pos) < distance) {
                tmp_aln->diff = split_events[i].pos - ref_pos;
                tmp_aln->high_error_side = false;
                if (split_events[i].strand && i > 0
                    && abs(split_events[i].read_pos_start - split_events[i-1].read_pos_stop) >= distance) {
                    tmp_aln->bp_read_pos = split_events[i].read_pos_start;
                    return true;
                } else if (!split_events[i].strand && i < split_events.size() - 1
                           && abs(split_events[i].read_pos_stop - split_events[i+1].read_pos_start) >= distance) {
                    tmp_aln->bp_read_pos = read_length - 1 - split_events[i].read_pos_stop;
                    return true;
                } else return false;
            } else if (abs(split_events[i].pos + split_events[i].length - ref_pos) < distance) {
                tmp_aln->diff = split_events[i].pos + split_events[i].length - ref_pos;
                tmp_aln->high_error_side = true;
                if (split_events[i].strand && i < split_events.size() - 1
                    && abs(split_events[i].read_pos_stop - split_events[i+1].read_pos_start) >= distance) {
                    tmp_aln->bp_read_pos = split_events[i].read_pos_stop;
                    return true;
                } else if (!split_events[i].strand && i > 0
                           && abs(split_events[i].read_pos_start - split_events[i-1].read_pos_stop) >= distance) {
                    tmp_aln->bp_read_pos =  read_length - 1 - split_events[i].read_pos_start;
                    return true;
                }
                else return false;
            }
        }
    }
    return false;
}

void detect_clipped_reads_rln(BpRln bp_realign,  Alignment* tmp_aln,
                              RefVector ref, std::map<std::string, ClippedRead> &mapClippedRead, ofstream& fasta_out){

    int distance = min(100, Parameter::Instance()->min_length);
    std::vector <aln_str> split_events = tmp_aln->getSA(ref);
    if (split_events.size() > Parameter::Instance()->max_splits)
        return;

    if (bp_realign.chr_pos.first - 20 <= tmp_aln->getPosition() ||
        bp_realign.chr_pos.first + 20 >= tmp_aln->getPosition() + tmp_aln->getRefLength()) {
        auto map_support = bp_realign.bp->get_coordinates().support;
        if (bp_realign.denovo){
            if (tmp_aln->getAlignment()->IsPrimaryAlignment() && !(tmp_aln->getAlignment()->AlignmentFlag & 0x800)) {
                if (map_support.find(tmp_aln->getName()) != map_support.end()) {
                    read_str read = map_support[tmp_aln->getName()];
//                clock_t begin_write = clock();
                    if (!read.processed) {
                        fasta_out << ">" << tmp_aln->getName() << "|" << read.coordinates.first
                                  << "|" << read.coordinates.second << endl;
                        fasta_out << tmp_aln->getQueryBases().substr(read.read_pos.first, read.read_pos.second -
                                                                                          read.read_pos.first) << endl;
                        read.processed = true;
                        bp_realign.bp->add_read(read, tmp_aln->getName());
                    }
//                Parameter::Instance()->meassure_time(begin_write, "time for writing reads");
                }
            }
            return;
        }

        if (map_support.find(tmp_aln->getName()) != map_support.end()) return;
        if (add_missed_tra(split_events, bp_realign, tmp_aln, ref)) return;
        
        bool existsGap = detect_gap(tmp_aln, split_events, bp_realign.chr_pos.first, distance,  ref);
        if (!existsGap) return;
        rln_str event_rln;
        event_rln.side = tmp_aln->high_error_side;
        event_rln.ref = bp_realign.chr;
        event_rln.ref_pos = bp_realign.chr_pos;
        event_rln.coordinate = bp_realign.coordinate;
        event_rln.bp = bp_realign.bp;
        event_rln.strand = bp_realign.isSameStrand;
//        if (tmp_aln->high_error_side)
        event_rln.read_pos = tmp_aln->bp_read_pos;
//        else event_rln.read_pos = tmp_aln->bp_read_pos - distance;
        event_rln.diff = tmp_aln->diff;

        if (mapClippedRead.find(tmp_aln->getName()) == mapClippedRead.end()){
            ClippedRead clipped_read;
            clipped_read.setChrId(tmp_aln->getRefID());
            clipped_read.set_ref_pos(split_events);
            if (tmp_aln->getAlignment()->IsPrimaryAlignment() && !(tmp_aln->getAlignment()->AlignmentFlag & 0x800) )
                clipped_read.setBases(tmp_aln->getQueryBases());
            clipped_read.events_rln.push_back(event_rln);
            mapClippedRead[tmp_aln->getName()] = clipped_read;
        } else {
            if (tmp_aln->getAlignment()->IsPrimaryAlignment() && !(tmp_aln->getAlignment()->AlignmentFlag & 0x800) )
                mapClippedRead[tmp_aln->getName()].setBases(tmp_aln->getQueryBases());
            mapClippedRead[tmp_aln->getName()].events_rln.push_back(event_rln);
        }
    }
}

void realign_across_read(BpRln bp_realign, vector<differences_str> &event_aln, Alignment* tmp_aln,
                         const bioio::FastaIndex  index, std::ifstream & fasta, RefVector ref, ofstream& fasta_out) {
    int distance = min(100, Parameter::Instance()->min_length);
    if (bp_realign.chr_pos.first - distance > tmp_aln->getPosition() &&
        bp_realign.chr_pos.first + distance < tmp_aln->getPosition() + tmp_aln->getRefLength()) {
        if (bp_realign.denovo ){
            if (tmp_aln->getAlignment()->IsPrimaryAlignment() && !(tmp_aln->getAlignment()->AlignmentFlag & 0x800)) {
                auto support_map = bp_realign.bp->get_coordinates().support;
                if (support_map.find(tmp_aln->getName()) != support_map.end()) {
                    read_str read = support_map[tmp_aln->getName()];
                    if (!read.processed) {
                        fasta_out << ">" << tmp_aln->getName() << "|" << read.coordinates.first
                                  << "|" << read.coordinates.second << endl;
                        fasta_out << tmp_aln->getQueryBases().substr(read.read_pos.first, read.read_pos.second -
                                                                                          read.read_pos.first) << endl;
                        read.processed = true;
                        bp_realign.bp->add_read(read, tmp_aln->getName());
                    }
                }
            }
            return;
        }

        bool exists_high_error_side = cal_high_error_side(event_aln, bp_realign.chr_pos.first, distance,  tmp_aln);
        if (!exists_high_error_side) return;
        int alt_aln_score = map_read(tmp_aln, bp_realign,  distance, index, fasta);
        if (alt_aln_score - tmp_aln->high_error_region_score > distance / 5 &&
            alt_aln_score > -0.2 * distance) {
            add_realign_read(bp_realign,  tmp_aln);
        }
    }
}

vector<aln_str> get_event(Alignment *& tmp){
    vector<aln_str> entries;
    aln_str event;
    event.RefID = tmp->getRefID();
    event.cigar = tmp->getCigar();
    event.length = (long) tmp->get_length(event.cigar);
    event.mq = tmp->getMappingQual();
    event.pos = (long) tmp->getPosition(); //+get_ref_lengths(event.RefID, ref);
    event.strand = tmp->getStrand();
    event.isMain = true;
    uint32_t sv;
    event.cross_N = ((sv & Ns_CLIPPED));

    tmp->get_coords(event, event.read_pos_start, event.read_pos_stop);

    entries.push_back(event);
   return entries;
}

void add_clipped_reads(Alignment *& tmp, std::vector<aln_str> events, short type, RefVector ref, IntervallTree& bst_rln, TNode *&root_rln, long read_id) {
    if (events.empty())
        events = get_event(tmp);
    int distance = min(100, Parameter::Instance()->min_length);
    long read_length = tmp->getAlignment()->Length;
    aln_str last_event = events[events.size()-1];

    if (events[0].read_pos_start >= distance) {
        position_str svs;
        read_str read;
        read.sequence = "NA";
        read.type = type;
        read.SV = 0;
        read.read_strand.first = events[0].strand;
        read.read_strand.second = events[0].strand;
        if (events[0].strand) {
            read.clipped_end = -1;
            read.coordinates.second = events[0].pos + get_ref_lengths(events[0].RefID, ref);
            read.coordinates.first = events[0].pos - distance + get_ref_lengths(events[0].RefID, ref);
        } else {
            read.clipped_end = 1;
            read.coordinates.first = events[0].pos + events[0].length + get_ref_lengths(events[0].RefID, ref);
            read.coordinates.second = events[0].pos + events[0].length + get_ref_lengths(events[0].RefID, ref)
                    + distance;
        }

        if (tmp->getStrand()){
//            read.read_pos.first = even?ts[0].read_pos_start - distance;
            read.read_pos.first = 0;
            read.read_pos.second =  events[0].read_pos_start;
        } else {
            read.read_pos.first = read_length - 1 - events[0].read_pos_start;
//            read.read_pos.second =  read_length - 1 - events[0].read_pos_start + distance;
            read.read_pos.second =  read_length - 1;
        }
//        cout << tmp->getName() << " " << tmp->getQueryBases().size() << " " << read.coordinates.first <<
//             " " << read.coordinates.second << endl;

        read.id = read_id;
        svs.start.min_pos = read.coordinates.first;
        svs.start.max_pos = read.coordinates.first;
        svs.stop.min_pos = read.coordinates.second;
        svs.stop.max_pos = read.coordinates.second;

        svs.support[tmp->getName()] = read;
        svs.support[tmp->getName()].length = abs(read.coordinates.second - read.coordinates.first);
        Breakpoint *point = new Breakpoint(svs, abs(read.coordinates.second - read.coordinates.first));
        bst_rln.insert(point, root_rln);
    }

    if (abs(last_event.read_pos_stop - read_length) > distance) {
        position_str svs;
        read_str read;
        read.sequence = "NA";
        read.type = type;
        read.SV = 0;
        read.read_strand.first = last_event.strand;
        read.read_strand.second = last_event.strand;
        if (!last_event.strand) {
            read.clipped_end = -1;
            read.coordinates.second = last_event.pos + get_ref_lengths(last_event.RefID, ref);
            read.coordinates.first = last_event.pos - distance+ get_ref_lengths(last_event.RefID, ref);
        } else {
            read.clipped_end = 1;
            read.coordinates.first = last_event.pos + last_event.length + get_ref_lengths(last_event.RefID, ref);
            read.coordinates.second = last_event.pos + last_event.length + distance + get_ref_lengths(last_event.RefID, ref);

        }
        if (tmp->getStrand()){
            read.read_pos.first = last_event.read_pos_stop;
//            read.read_pos.second =  last_event.read_pos_stop + distance;
            read.read_pos.second = read_length - 1;
        } else {
//            read.read_pos.first = read_length - 1 - last_event.read_pos_stop - distance;
            read.read_pos.first = 0;
            read.read_pos.second = read_length - 1 - last_event.read_pos_stop;
        }
//        cout << tmp->getName() << " " << tmp->getQueryBases().size() << " " << read.coordinates.first <<
//                " " << read.coordinates.second << endl;

        read.id = read_id;
        svs.start.min_pos = read.coordinates.first;
        svs.start.max_pos = read.coordinates.first;
        svs.stop.min_pos = read.coordinates.second;
        svs.stop.max_pos = read.coordinates.second;

        svs.support[tmp->getName()] = read;
        svs.support[tmp->getName()].length = abs(read.coordinates.second - read.coordinates.first);
        Breakpoint *point = new Breakpoint(svs, abs(read.coordinates.second - read.coordinates.first));
        bst_rln.insert(point, root_rln);
    }


    for (size_t i = 1; i < events.size(); i++) {
        //	if (events[i - 1].length > Parameter::Instance()->min_segment_size || events[i].length > Parameter::Instance()->min_segment_size) {
        position_str svs;
        //position_str stop;
        read_str read;
        read.sequence = "NA";
        //read.name = tmp->getName();
        read.type = type;
        read.SV = 0;
        read.read_strand.first = events[i - 1].strand;
        read.read_strand.second = !events[i].strand;
        if (events[i-1].RefID == events[i].RefID) {
            if (events[i].read_pos_start - events[i - 1].read_pos_stop >= distance) {
                if (events[i - 1].strand == events[i].strand) {
                    //check this with + - strands!!

                    if (events[i - 1].strand) { //"++"
                        svs.start.min_pos = events[i - 1].pos + events[i - 1].length + get_ref_lengths(events[i - 1].RefID, ref);
                        svs.stop.max_pos = events[i].pos + get_ref_lengths(events[i].RefID, ref);
                    } else { //"--"
                        svs.start.min_pos = events[i - 1].pos + get_ref_lengths(events[i - 1].RefID, ref);
                        svs.stop.max_pos = events[i].pos + events[i].length + get_ref_lengths(events[i].RefID, ref);
                    }

                } else {
                    if (events[i - 1].strand) { //"+-"
                        svs.start.min_pos =
                                events[i - 1].pos + events[i - 1].length + get_ref_lengths(events[i - 1].RefID, ref);
                        svs.stop.max_pos = events[i].pos + events[i].length + get_ref_lengths(events[i].RefID, ref);
                    } else { //"-+"
                        svs.start.min_pos = events[i - 1].pos + get_ref_lengths(events[i - 1].RefID, ref);
                        svs.stop.max_pos = events[i].pos + get_ref_lengths(events[i].RefID, ref);
                    }
                }
                if (tmp->getStrand()) {
                    read.read_pos.first = events[i - 1].read_pos_stop;
                    read.read_pos.second = events[i].read_pos_start;
                } else {
                    read.read_pos.first = read_length - 1 - events[i].read_pos_start;
                    read.read_pos.second = read_length - 1 - events[i - 1].read_pos_stop;
                }

                read.SV |= TRA;



                //std::cout<<"split"<<std::endl;

                svs.start.max_pos = svs.start.min_pos;
                svs.stop.min_pos = svs.stop.max_pos;
                if (svs.start.min_pos > svs.stop.max_pos) {
                    //maybe we have to invert the directions???
                    svs_breakpoint_str pos = svs.start;
                    svs.start = svs.stop;
                    svs.stop = pos;

                    pair<bool, bool> tmp = read.strand;

                    read.strand.first = tmp.second;
                    read.strand.second = tmp.first;
                }

                //TODO: we might not need this:
                if (svs.start.min_pos > svs.stop.max_pos) {
                    read.coordinates.first = svs.stop.max_pos;
                    read.coordinates.second = svs.start.min_pos;
                } else {
                    read.coordinates.first = svs.start.min_pos;
                    read.coordinates.second = svs.stop.max_pos;
                }
//                cout << tmp->getName() << " " << tmp->getQueryBases().size() << " " << read.coordinates.first <<
//                     " " << read.coordinates.second << endl;

                //pool out?
                read.id = read_id;
                svs.support[tmp->getName()] = read;
                svs.support[tmp->getName()].length = abs(read.coordinates.second - read.coordinates.first);
                Breakpoint *point = new Breakpoint(svs, abs(read.coordinates.second - read.coordinates.first));
                //std::cout<<"split ADD: " << <<" Name: "<<tmp->getName()<<" "<< svs.start.min_pos- get_ref_lengths(events[i].RefID, ref)<<"-"<<svs.stop.max_pos- get_ref_lengths(events[i].RefID, ref)<<std::endl;
                bst_rln.insert(point, root_rln);

            }
        }
    }
}

bool add_missed_tra(vector<aln_str> split_events, BpRln bpRln, Alignment * tmp_aln, RefVector ref){
    int tol = 20;
    for (int i = 1; i < split_events.size(); i++) {
        read_str read;
        if (split_events[i-1].RefID == bpRln.chr_idx.first && abs(split_events[i-1].pos - bpRln.chr_pos.first) <= tol
            && split_events[i].RefID == bpRln.chr_idx.second && abs(split_events[i].pos - bpRln.chr_pos.second) <= tol
            && split_events[i-1].mq >= Parameter::Instance()->min_mq && split_events[i].mq >= Parameter::Instance()->min_mq) {
            read.coordinates.first = split_events[i-1].pos + get_ref_lengths(split_events[i-1].RefID, ref);
            read.coordinates.second = split_events[i].pos + get_ref_lengths(split_events[i].RefID, ref);
        } else if (split_events[i].RefID == bpRln.chr_idx.first && abs(split_events[i].pos - bpRln.chr_pos.first) <= tol
            && split_events[i-1].RefID == bpRln.chr_idx.second && abs(split_events[i-1].pos - bpRln.chr_pos.second) <= tol
            && split_events[i-1].mq >= Parameter::Instance()->min_mq && split_events[i].mq >= Parameter::Instance()->min_mq) {
            read.coordinates.first = split_events[i].pos + get_ref_lengths(split_events[i].RefID, ref);
            read.coordinates.second = split_events[i-1].pos + get_ref_lengths(split_events[i-1].RefID, ref);

        } else continue;
        if (read.coordinates.first > read.coordinates.second){
            read_str tmp = read;
            read.coordinates.first = tmp.coordinates.second;
            read.coordinates.second = tmp.coordinates.first;
        }
        read.SV = TRA;
        read.type = 1;
        bpRln.bp->add_read(read, tmp_aln->getName());
        return true;
    }
    return false;
}

void add_high_error_reads(Alignment *& tmp_aln, RefVector ref, IntervallTree& bst_rln, TNode *&root_rln, long read_id){

    std::vector<indel_str> dels;
    vector<differences_str> event_aln;

    event_aln = tmp_aln->summarizeAlignment(dels);
    for (size_t j = 0; j < event_aln.size(); j++) {
        if (event_aln[j].type == 0)
            fix_mismatch_read_pos(event_aln, j, tmp_aln);
    }

    PlaneSweep_slim * plane = new PlaneSweep_slim();
    vector<pair_str> profile;
    for (size_t i = 0; i < event_aln.size(); i++) {
        pair_str tmp;
        tmp.position = -1;
        if (event_aln[i].type == 0) { //subst_rlnitutions.
            tmp = plane->add_mut(event_aln[i].position, 1, Parameter::Instance()->window_thresh);
        } else if (abs(event_aln[i].type) <= 10){
            tmp = plane->add_mut(event_aln[i].position, abs(event_aln[i].type), Parameter::Instance()->window_thresh);	// abs(event_aln[i].type)
        }
        if (tmp.position != -1 && (profile.empty() || (tmp.position - profile[profile.size() - 1].position) > 100)) {	//for noisy events;
            profile.push_back(tmp);
        }
    }
    if (profile.size() > 4)
        return;
    int stop = 0;
    size_t start = 0;
    for (size_t i = 0; i < profile.size() && stop < event_aln.size(); i++) {
        if (profile[i].position >= event_aln[stop].position) {
            //find the position:
            size_t pos = 0;
            while (pos < event_aln.size() && event_aln[pos].position != profile[i].position) {
                pos++;
            }
            //run back to find the start:
            start = pos;
            int prev = event_aln[pos].position;
            start = pos;
            int prev_type = 1;
            int start_pos = event_aln[pos].position;
            //todo it is actually pos + type and not *type
            while (start > 0 ) {    //13		//} * abs(event_aln[start].type) + 1)) { //TODO I  dont like 13!??
                if (prev - start_pos >= Parameter::Instance()->max_dist_alns) {
                    int leftSideErrors = 0;
                    for (int i = start+1; event_aln[start+1].position - event_aln[i].position < 50; i--){
                        leftSideErrors += max(1, (int)abs(event_aln[i].type));
                        if (i== 0) break;
                    }
                    if (leftSideErrors < 15) break;
                }
                prev = event_aln[start].position;
                start--;
                start_pos = event_aln[start].position + max(1, (int)abs(event_aln[start].type));
            }
            if (start + 1 < event_aln.size()) { //TODO do some testing!
                start++; //we are running one too far!
            }
            //run forward to identify the stop:
            prev = event_aln[pos].position;
            stop = pos;
            prev_type = 1;
            while (stop < event_aln.size()) {        // * abs(event_aln[stop].type) + 1)) {
                if (event_aln[stop].position - prev >= Parameter::Instance()->max_dist_alns){
                    int rightSideErrors = 0;
                    for (int i = stop-1; event_aln[i].position + abs(event_aln[i].type) - event_aln[stop-1].position < 50; i++){
                        rightSideErrors += max(1, (int)abs(event_aln[i].type));
                        if (i >= event_aln.size() - 1) break;
                    }
                    if (rightSideErrors < 15) break;
                }
                prev = event_aln[stop].position;
                prev_type = abs(event_aln[stop].type);
                if (prev_type == 0) {
                    prev_type = 1;
                }
                prev += prev_type;
                stop++;
            }
            if (stop > 0) {
                stop--;
            }
//            cout << event_aln[start].position << " " << event_aln[stop].position << endl;
            position_str svs;
            read_str read;
            read.sequence = "NA";
            read.type = 1;
            read.SV = TRA;
            read.read_strand.first = tmp_aln->getStrand();
            read.read_strand.second =  tmp_aln->getStrand();
            read.coordinates.first = event_aln[start].position + get_ref_lengths(tmp_aln->getRefID(), ref);
            read.coordinates.second = event_aln[stop].position + get_ref_lengths(tmp_aln->getRefID(), ref);
            read.read_pos.first = event_aln[start].readposition;
            read.read_pos.second = event_aln[stop].readposition;
            read.id = read_id;
            svs.start.min_pos = read.coordinates.first;
            svs.start.max_pos = read.coordinates.first;
            svs.stop.min_pos = read.coordinates.second;
            svs.stop.max_pos = read.coordinates.second;
            svs.support[tmp_aln->getName()] = read;
            svs.support[tmp_aln->getName()].length = abs(read.coordinates.second - read.coordinates.first);
            Breakpoint *point = new Breakpoint(svs, abs(read.coordinates.second - read.coordinates.first));
            bst_rln.insert(point, root_rln);
//            cout << ">" << tmp_aln->getName() << endl;
//            cout << tmp_aln->getQueryBases().substr(read.read_pos.first, read.read_pos.second - read.read_pos.first)
//                 << endl;
        }
    }
}


void detect_bps_for_realn(vector<Breakpoint*> points, vector<Breakpoint*> points_rln, const RefVector ref, vector<BpRln> &bp_rlns) {
    long max_dist = Parameter::Instance()->max_dist * 2;
    vector<Breakpoint*> points_rln_final,points_final;
    for (Breakpoint* bp_rln: points_rln){
        if (bp_rln->get_coordinates().support.size() >= max(1, Parameter::Instance()->min_support / 2))
//            if (bp_rln->get_coordinates().support.size() >= 1)
            points_rln_final.push_back(bp_rln);
        else {
            Breakpoint * a = bp_rln;
            delete a;
            a = NULL;
        }
    }
    points_rln.clear();

    if (points.empty()){
        for (auto bp: points_rln_final) detect_bp_for_realn(bp, ref, bp_rlns, true);
        return;
    }

    for (auto bp: points){
        if (bp->get_SVtype() & TRA)
            detect_bp_for_realn(bp, ref, bp_rlns, false);
        else if (bp->get_coordinates().support.size() < Parameter::Instance()->min_support)
             points_final.push_back(bp);
    }
    points.clear();
    for (Breakpoint* bp_rln: points_rln_final) bp_rln->set_valid(true);
//    Breakpoint * bp;
//    bp->get_coordinates().start.min_pos

    size_t j = 0;

    for (size_t i = 0; i < points_final.size(); i++) {

        if (points_final[i]->get_coordinates().start.min_pos <
            points_rln_final[j]->get_coordinates().start.min_pos - max_dist)
            continue;

        while (points_final[i]->get_coordinates().start.min_pos >
               points_rln_final[j]->get_coordinates().start.min_pos + max_dist && j < points_rln_final.size()) {
            j++;
        }
        if (j == points_rln_final.size())
            break;

        for (size_t k = j; points_rln_final[k]->get_coordinates().start.min_pos <= points_final[i]->get_coordinates().start.min_pos + max_dist
        && k < points_rln_final.size(); k++) {
            if (!points_rln_final[k]->get_valid()) continue;
            vector<string> same_reads;
            set<string> i_reads, k_reads;
            auto i_support_map = points_final[i]->get_coordinates().support;
            auto k_support_map =  points_rln_final[k]->get_coordinates().support;
            for (auto it = i_support_map.begin(); it != i_support_map.end(); it++) i_reads.insert(it->first);
            for (auto it = k_support_map.begin(); it != k_support_map.end(); it++) k_reads.insert(it->first);
            set_intersection(i_reads.begin(),i_reads.end(),k_reads.begin(),k_reads.end(), std::back_inserter(same_reads));

            int num_overlap = 0;
            for (string read_name: same_reads){
                read_str i_read = points_final[i]->get_coordinates().support[read_name];
                read_str k_read = points_final[k]->get_coordinates().support[read_name];

                 if ( i_read.coordinates.first <= k_read.coordinates.second &&
                    i_read.coordinates.second >= k_read.coordinates.first )
                    num_overlap++;
            }
            if (num_overlap > max(1, Parameter::Instance()->min_support / 2))
                points_rln_final[j]->set_valid(false);
        }
    }

    std::sort(bp_rlns.begin(), bp_rlns.end());
    vector<BpRln> bps_tmp;

    for (Breakpoint * bp: points_rln_final) {
        if (bp->get_valid()) {
            detect_bp_for_realn(bp, ref, bps_tmp, true);
        }
    }
    size_t point_size = bp_rlns.size(), l = 0;

    for (int i = 0; i < point_size; i++) {
        if (bps_tmp.size() == 0)
            break;
        while (bp_rlns[i].coordinate.first > bps_tmp[l].coordinate.second + Parameter::Instance()->min_length &&
        l < bps_tmp.size()) {
            bp_rlns.push_back(bps_tmp[l]);
            pair<long, long> coordinate = bps_tmp[l].coordinate;
            std::swap(coordinate.first, coordinate.second);
            BpRln stop_bp(bps_tmp[l].isSameStrand, coordinate, ref, bps_tmp[l].bp);
            stop_bp.denovo = true;
            bp_rlns.push_back(stop_bp);
            l++;
        }
        while ((bp_rlns[i].coordinate.first >= bps_tmp[l].coordinate.first - Parameter::Instance()->min_length &&
                bp_rlns[i].coordinate.first <= bps_tmp[l].coordinate.second + Parameter::Instance()->min_length) &&
               l < bps_tmp.size())
            l++;
        if (l == bps_tmp.size())
            break;
    }

    for (; l < bps_tmp.size(); l++){
        bp_rlns.push_back(bps_tmp[l]);
        pair<long, long> coordinate = bps_tmp[l].coordinate;
        std::swap(coordinate.first, coordinate.second);
        BpRln stop_bp(bps_tmp[l].isSameStrand, coordinate, ref, bps_tmp[l].bp);
        stop_bp.denovo = true;
        bp_rlns.push_back(stop_bp);
        }
    std::sort(bp_rlns.begin(), bp_rlns.end());
}

int getRefIdx(RefVector ref, string chr){
    for (int i = 0; i < ref.size(); i++){
        if (ref[i].RefName == chr)
            return i;
    }
}

Breakpoint * generate_bp(pair<long, long> coordinate, string name) {
    read_str read;
    read.sequence = "NA";
    read.type = 1;
    read.SV = TRA;
    read.read_strand.first = true;
    read.read_strand.second = true;
    position_str svs;
    svs.start.min_pos = coordinate.first;
    svs.stop.max_pos = coordinate.second;
    svs.start.max_pos = svs.start.min_pos;
    svs.stop.min_pos = svs.stop.max_pos;
    if (svs.start.min_pos > svs.stop.max_pos) {
        //maybe we have to invert the directions???
        svs_breakpoint_str pos = svs.start;
        svs.start = svs.stop;
        svs.stop = pos;

        pair<bool, bool> tmp = read.strand;

        read.strand.first = tmp.second;
        read.strand.second = tmp.first;
    }

    //TODO: we might not need this:
    if (svs.start.min_pos > svs.stop.max_pos) {
        read.coordinates.first = svs.stop.max_pos;
        read.coordinates.second = svs.start.min_pos;
    } else {
        read.coordinates.first = svs.start.min_pos;
        read.coordinates.second = svs.stop.max_pos;
    }
    svs.support[name] = read;
    svs.support[name].length = abs(read.coordinates.second - read.coordinates.first);
    Breakpoint * point = new Breakpoint(svs, abs(read.coordinates.second - read.coordinates.first));
    return point;

}

void minimap2(string target_path, string query_path, RefVector ref, IntervallTree &bst_rln_rln_denovo, TNode * &root_rln_rln_denovo){
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    int n_threads = 10;

    int distance = min(100, Parameter::Instance()->min_length);
    mm_verbose = 2; // disable message output to stderr
    mm_set_opt(0, &iopt, &mopt);
    mopt.flag |= MM_F_CIGAR; // perform alignment

    // open query file for reading; you may use your favorite FASTA/Q parser
    gzFile f = gzopen(query_path.c_str(), "r");
    assert(f);
    kseq_t *ks = kseq_init(f);

    // open index reader
    mm_idx_reader_t *r = mm_idx_reader_open(target_path.c_str(), &iopt, 0);
    mm_idx_t *mi;
    cout << "start global remapping..." << endl;
    while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
        mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!
        mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
        while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
            mm_reg1_t *reg;
            int j, i, n_reg;
            reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query
            for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
                mm_reg1_t *r = &reg[j];
                assert(r->p); // with MM_F_CIGAR, this should not be NULL
                string str(ks->name.s);
                string chr(mi->seq[r->rid].name);
                vector<string> items;
                split(str, '|', items);
                string name = items[0];
                pair<long, long> coordinate;
                int refIdx = getRefIdx(ref, chr);
                coordinate.first = std::stol(items[1]);
                coordinate.second = std::stol(items[2]);
                if (coordinate.first > coordinate.second)
                    std::swap(coordinate.first, coordinate.second);
                if ("+-"[r->rev] == '-')
                    std::swap(coordinate.first, coordinate.second);
                int query_length = ks->seq.l, map_quality = r->mapq;
                if (map_quality > 3 && r->qs < 20 && r->qe > query_length - 20 && r->mlen > r->blen * 0.8){
                    if (abs(coordinate.second - coordinate.first) > distance){
                        position_str svs;
                        //position_str stop;
                        pair<long, long> coordinate_tra;
                        coordinate_tra.first = coordinate.first;
                        coordinate_tra.second = r->rs + get_ref_lengths(refIdx, ref);
                        Breakpoint* point0 = generate_bp(coordinate_tra, name);
                        bst_rln_rln_denovo.insert(point0, root_rln_rln_denovo);
                        coordinate_tra.first = coordinate.second;
                        coordinate_tra.second = r->re + get_ref_lengths(refIdx, ref);
                        Breakpoint* point1 = generate_bp(coordinate_tra, name);
                        bst_rln_rln_denovo.insert(point1, root_rln_rln_denovo);
                    }
                }
//                printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
//                printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:", mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq);
//                for (i = 0; i < r->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
//                    printf("%d%c", r->p->cigar[i]>>4, "MIDNSH"[r->p->cigar[i]&0xf]);
                putchar('\n');
                free(r->p);
            }
            free(reg);
        }
        mm_tbuf_destroy(tbuf);
        mm_idx_destroy(mi);
    }
    cout << "finish global remapping" << endl;
    mm_idx_reader_close(r); // close the index reader
    kseq_destroy(ks); // close the query file
    gzclose(f);
}

