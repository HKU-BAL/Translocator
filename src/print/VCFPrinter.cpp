/*
 * VCFPrinter.cpp
 *
 *  Created on: Aug 24, 2015
 *      Author: fsedlaze
 */

#include "VCFPrinter.h"

void VCFPrinter::print_header() {
	fprintf(file, "%s", "##fileformat=VCFv4.1\n");
	fprintf(file, "%s", "##source=Sniffles\n");
	string time = currentDateTime();
	fprintf(file, "%s", "##fileDate=");
	fprintf(file, "%s", time.c_str());

	//REport over all chrs:
	for (size_t i = 0; i < this->ref.size(); i++) {
		fprintf(file, "%s", "\n");
		fprintf(file, "%s", "##contig=<ID=");
		fprintf(file, "%s", ref[i].RefName.c_str());
		fprintf(file, "%s", ",length=");
		fprintf(file, "%i", (int) ref[i].RefLength);
		fprintf(file, "%c", '>');
	}

	fprintf(file, "%s", "\n");
	fprintf(file, "%s", "##ALT=<ID=DEL,Description=\"Deletion\">\n");
	fprintf(file, "%s", "##ALT=<ID=DUP,Description=\"Duplication\">\n");
	fprintf(file, "%s", "##ALT=<ID=INV,Description=\"Inversion\">\n");
	fprintf(file, "%s", "##ALT=<ID=INVDUP,Description=\"InvertedDUP with unknown boundaries\">\n");
	fprintf(file, "%s", "##ALT=<ID=TRA,Description=\"Translocation\">\n");
	fprintf(file, "%s", "##ALT=<ID=INS,Description=\"Insertion\">\n");
	fprintf(file, "%s", "##FILTER=<ID=UNRESOLVED,Description=\"An insertion that is longer than the read and thus we cannot predict the full size.\">\n");
	fprintf(file, "%s", "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n");
	fprintf(file, "%s", "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n");
	fprintf(file, "%s", "##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">\n");
	fprintf(file, "%s", "##INFO=<ID=RE,Number=1,Type=Integer,Description=\"read support\">\n");
	fprintf(file, "%s", "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n");
	fprintf(file, "%s", "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n");

	//##FILTER=<ID=LowQual,Description="PE/SR support below 3 or mapping quality below 20.">
	fprintf(file, "%s", "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n");
	fprintf(file, "%s", "##INFO=<ID=REF_strand,Number=2,Type=Integer,Description=\"Length of the SV\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">\n");
	fprintf(file, "%s", "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
	if (Parameter::Instance()->report_n_reads > 0 || Parameter::Instance()->report_n_reads == -1) {
		fprintf(file, "%s", "##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"Names of reads supporting SVs (comma separated)\">\n");
	}
	if (Parameter::Instance()->print_seq) {
		fprintf(file, "%s", "##INFO=<ID=SEQ,Number=1,Type=String,Description=\"Extracted sequence from the best representative read.\">\n");
	}

	if (Parameter::Instance()->read_strand) {
		fprintf(file, "%s", "##INFO=<ID=STRANDS2,Number=4,Type=Integer,Description=\"alt reads first + ,alt reads first -,alt reads second + ,alt reads second -.\">\n");
		fprintf(file, "%s", "##INFO=<ID=REF_strand,Number=2,Type=Integer,Description=\"plus strand ref, minus strand ref.\">\n");
	}

	fprintf(file, "%s", "##INFO=<ID=STD_quant_start,Number=A,Type=Float,Description=\"STD of the start breakpoints across the reads.\">\n");
	fprintf(file, "%s", "##INFO=<ID=STD_quant_stop,Number=A,Type=Float,Description=\"STD of the stop breakpoints across the reads.\">\n");
	fprintf(file, "%s", "##INFO=<ID=Kurtosis_quant_start,Number=A,Type=Float,Description=\"Kurtosis value of the start breakpoints across the reads.\">\n");
	fprintf(file, "%s", "##INFO=<ID=Kurtosis_quant_stop,Number=A,Type=Float,Description=\"Kurtosis value of the stop breakpoints across the reads.\">\n");
	fprintf(file, "%s", "##INFO=<ID=SUPTYPE,Number=A,Type=String,Description=\"Type by which the variant is supported.(SR,ALN,NR)\">\n");
	fprintf(file, "%s", "##INFO=<ID=STRANDS,Number=A,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)\">\n");
	fprintf(file, "%s", "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency.\">\n");
	fprintf(file, "%s", "##INFO=<ID=ZMW,Number=A,Type=Integer,Description=\"Number of ZMWs (Pacbio) supporting SV.\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference reads\">\n");
	fprintf(file, "%s", "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant reads\">\n");

	fprintf(file, "%s", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
	for (size_t i = 0; i < Parameter::Instance()->bam_files.size(); i++) {
		fprintf(file, "%c", '\t');
		fprintf(file, "%s", Parameter::Instance()->bam_files[i].c_str());
	}
	fprintf(file, "%c", '\n');

}
void VCFPrinter::print_body(Breakpoint * &SV, RefVector ref) {
	if (!this->bed_tree.is_in(SV->get_coordinates().start.most_support, this->root) && !this->bed_tree.is_in(SV->get_coordinates().stop.most_support, this->root)) {
		//temp. store read names supporting this SVs to later group the SVs together.
		double std_quant_start = 0;
		double std_quant_stop = 0;

		pair<double, double> kurtosis;
		pair<double, double> std_quant;
		double std_length = 0;
		int zmws = 0;
		bool ok_to_print = (to_print(SV, std_quant, kurtosis, std_length, zmws) || Parameter::Instance()->ignore_std);
		//std::cout << "Print check: " << std_quant.first << " " << std_quant.second << endl;
		if (ok_to_print && (zmws == 0 || zmws >= Parameter::Instance()->min_zmw)) {
			if (Parameter::Instance()->phase) {
				store_readnames(SV->get_read_ids(), id);
			}
			std::string chr;
			int start = IPrinter::calc_pos(SV->get_coordinates().start.most_support, ref, chr);
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", '\t');
			if (start < 1) {
				start = 1;
			}
			fprintf(file, "%i", start);
			fprintf(file, "%c", '\t');
			fprintf(file, "%i", id);
			id++;

			long end_coord = SV->get_coordinates().stop.most_support;
			if (((SV->get_SVtype() & INS))) { // && SV->get_length() == Parameter::Instance()->huge_ins) && SV->get_types().is_ALN) {
				end_coord = std::max((SV->get_coordinates().stop.most_support - (long) SV->get_length()), (long) start);
			}

			int end = IPrinter::calc_pos(end_coord, ref, chr);
			if (end < 1) {
				end = 1;
			}
			std::string strands = SV->get_strand(1);

			if (Parameter::Instance()->reportBND && (SV->get_SVtype() & TRA)) {
				//N[22:36765684[ +-
				//]21:10540232]N -+
				fprintf(file, "%s", "\tN\t");
				if (strands[0] == '-') { //&&
					fprintf(file, "%s", "]");
					fprintf(file, "%s", chr.c_str());
					fprintf(file, "%c", ':');
					fprintf(file, "%i", end);
					fprintf(file, "%s", "]N");

				} else {
					fprintf(file, "%s", "N[");
					fprintf(file, "%s", chr.c_str());
					fprintf(file, "%c", ':');
					fprintf(file, "%i", end);
					fprintf(file, "%c", '[');
				}
			} else if (!SV->get_sequence().empty() && ((SV->get_SVtype() & INS) || (SV->get_SVtype() & DEL))) {
				fprintf(file, "%c", '\t');
				if ((SV->get_SVtype() & DEL)) {
					fprintf(file, "%s", SV->get_sequence().c_str());
				} else {
					fprintf(file, "%c", 'N');
				}
				fprintf(file, "%c", '\t');
				if ((SV->get_SVtype() & INS)) {
					fprintf(file, "%s", SV->get_sequence().c_str());
				} else {
					fprintf(file, "%c", 'N');
				}

			} else {
				fprintf(file, "%s", "\tN\t");
				fprintf(file, "%c", '<');
				fprintf(file, "%s", IPrinter::get_type(SV->get_SVtype()).c_str());
				fprintf(file, "%c", '>');
			}

			if (((SV->get_SVtype() & INS) && SV->get_length() == Parameter::Instance()->huge_ins) && SV->get_types().is_ALN) {
				fprintf(file, "%s", "\t.\tUNRESOLVED\t");
			} else {
				fprintf(file, "%s", "\t.\tPASS\t");
			}

			if (std_quant.first < 10 && std_quant.second < 10) {
				fprintf(file, "%s", "PRECISE");
			} else {
				fprintf(file, "%s", "IMPRECISE");
			}

			fprintf(file, "%s", ";SVMETHOD=Snifflesv");
			fprintf(file, "%s", Parameter::Instance()->version.c_str());

			if (!(Parameter::Instance()->reportBND && (SV->get_SVtype() & TRA))) {

				fprintf(file, "%s", ";CHR2=");
				fprintf(file, "%s", chr.c_str());
				fprintf(file, "%s", ";END=");

				if (SV->get_SVtype() & INS) {
					fprintf(file, "%i", std::max((int) end, start));
				} else {

					fprintf(file, "%i", end);
				}
			}
			if (zmws != 0) {
				fprintf(file, "%s", ";ZMW=");
				fprintf(file, "%i", zmws);
			}
			fprintf(file, "%s", ";STD_quant_start=");
			fprintf(file, "%f", std_quant.first);
			fprintf(file, "%s", ";STD_quant_stop=");
			fprintf(file, "%f", std_quant.second);
			fprintf(file, "%s", ";Kurtosis_quant_start=");
			fprintf(file, "%f", kurtosis.first);
			fprintf(file, "%s", ";Kurtosis_quant_stop=");
			fprintf(file, "%f", kurtosis.second);

			fprintf(file, "%s", ";SVTYPE=");
			if (Parameter::Instance()->reportBND && (SV->get_SVtype() & TRA)) {
				fprintf(file, "%s", "BND");
			} else {
				fprintf(file, "%s", IPrinter::get_type(SV->get_SVtype()).c_str());
			}

			if (Parameter::Instance()->report_n_reads > 0 || Parameter::Instance()->report_n_reads == -1) {
				fprintf(file, "%s", ";RNAMES=");
				fprintf(file, "%s", SV->get_read_names().c_str());
			}
			fprintf(file, "%s", ";SUPTYPE=");
			fprintf(file, "%s", SV->get_supporting_types().c_str());
			fprintf(file, "%s", ";SVLEN=");

			if (((SV->get_SVtype() & INS) && SV->get_length() == Parameter::Instance()->huge_ins) && SV->get_types().is_ALN) {
				if (SV->get_sequence().size() != 0) { //!
					fprintf(file, "%i", SV->get_sequence().size());
				} else {
					fprintf(file, "%i", 1);
				}
			} else if (SV->get_SVtype() & TRA) {
				fprintf(file, "%i", 1);
			} else if (SV->get_SVtype() & DEL) {
				fprintf(file, "%i", SV->get_length() * -1);
			} else {
				fprintf(file, "%i", SV->get_length());
			}
			//	}
			fprintf(file, "%s", ";STRANDS=");
			fprintf(file, "%s", strands.c_str());
			if (Parameter::Instance()->read_strand) {
				fprintf(file, "%s", ";STRANDS2=");
				std::map<std::string, read_str> support = SV->get_coordinates().support;
				pair<int, int> tmp_start;
				pair<int, int> tmp_stop;
				tmp_start.first = 0;
				tmp_start.second = 0;
				tmp_stop.first = 0;
				tmp_stop.second = 0;
				for (std::map<std::string, read_str>::iterator i = support.begin(); i != support.end(); i++) {
					if ((*i).second.read_strand.first) {
						tmp_start.first++;
					} else {
						tmp_start.second++;
					}
					if ((*i).second.read_strand.second) {
						tmp_stop.first++;
					} else {
						tmp_stop.second++;
					}
				}
				fprintf(file, "%i", tmp_start.first);
				fprintf(file, "%s", ",");
				fprintf(file, "%i", tmp_start.second);
				fprintf(file, "%s", ",");
				fprintf(file, "%i", tmp_stop.first);
				fprintf(file, "%s", ",");
				fprintf(file, "%i", tmp_stop.second);
			}

			//	if (Parameter::Instance()->print_seq && !SV->get_sequence().empty()) {
			//		fprintf(file, "%s", ";SEQ=");
			//		fprintf(file, "%s", SV->get_sequence().c_str());
			//	}
			fprintf(file, "%s", ";RE=");
			fprintf(file, "%i", SV->get_support());
			//if(Parameter::Instance()->genotype){
			fprintf(file, "%s", "\tGT:DR:DV\t./.:.:");
			fprintf(file, "%i", SV->get_support());
			//}else{
			//	fprintf(file, "%s",this->assess_genotype(SV->get_refcount(),SV->get_support()).c_str());
			//}
			fprintf(file, "%c", '\n');
		}
	}
}

void VCFPrinter::print_body_recall(Breakpoint * &SV, RefVector ref) {
	if (Parameter::Instance()->phase) {
		store_readnames(SV->get_read_ids(), id);
	}

	pair<double, double> kurtosis;
	pair<double, double> std_quant;
	double std_length = 0;
	int zmws = 0;
	bool ok_to_print = to_print(SV, std_quant, kurtosis, std_length, zmws);

	std::string chr;
	int start = IPrinter::calc_pos(SV->get_coordinates().start.most_support, ref, chr);
	fprintf(file, "%s", chr.c_str());
	fprintf(file, "%c", '\t');
	fprintf(file, "%i", start);
	fprintf(file, "%c", '\t');
	fprintf(file, "%i", id);
	id++;

	int end = IPrinter::calc_pos(SV->get_coordinates().stop.most_support, ref, chr);
	std::string strands = SV->get_strand(1);
	fprintf(file, "%s", "\tN\t");
	if (Parameter::Instance()->reportBND && (SV->get_SVtype() & TRA)) {
		//N[22:36765684[ +-
		//]21:10540232]N -+
		if (strands[0] == '-' && strands[0] == '+') {
			fprintf(file, "%s", "]");
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", ':');
			fprintf(file, "%i", end);
			fprintf(file, "%s", "]N");

		} else {
			fprintf(file, "%s", "N[");
			fprintf(file, "%s", chr.c_str());
			fprintf(file, "%c", ':');
			fprintf(file, "%i", end);
			fprintf(file, "%c", '[');
		}
	} else {

		fprintf(file, "%c", '<');
		fprintf(file, "%s", IPrinter::get_type(SV->get_SVtype()).c_str());
		fprintf(file, "%c", '>');
	}

	fprintf(file, "%s", "\t.\tPASS\t");
	fprintf(file, "%s", "IMPRECISE");
	fprintf(file, "%s", ";SVMETHOD=Snifflesv");
	fprintf(file, "%s", Parameter::Instance()->version.c_str());
	if (!(Parameter::Instance()->reportBND && (SV->get_SVtype() & TRA))) {
		fprintf(file, "%s", ";CHR2=");
		fprintf(file, "%s", chr.c_str());
		fprintf(file, "%s", ";END=");

		if (SV->get_SVtype() & INS) {
			fprintf(file, "%i", std::max((int) (end - SV->get_length()), start));
		} else {
			fprintf(file, "%i", end);
		}
	}

	fprintf(file, "%s", ";SVTYPE=");
	if (Parameter::Instance()->reportBND && (SV->get_SVtype() & TRA)) {
		fprintf(file, "%s", "BND");
	} else {
		fprintf(file, "%s", IPrinter::get_type(SV->get_SVtype()).c_str());
	}
	if (Parameter::Instance()->report_n_reads > 0 || Parameter::Instance()->report_n_reads == -1) {
		fprintf(file, "%s", ";RNAMES=");
		fprintf(file, "%s", SV->get_read_names().c_str());
	}
	fprintf(file, "%s", ";SUPTYPE=");
	fprintf(file, "%s", SV->get_supporting_types().c_str());
	fprintf(file, "%s", ";SVLEN=");

	if (((SV->get_SVtype() & INS) && SV->get_length() == Parameter::Instance()->huge_ins) && !SV->get_types().is_SR) {
		fprintf(file, "%s", "NA");
	} else {
		fprintf(file, "%i", SV->get_length());
	}
	//	}
	fprintf(file, "%s", ";STRANDS=");
	fprintf(file, "%s", strands.c_str());
	fprintf(file, "%s", ";SEQ=");
	fprintf(file, "%s", SV->get_sequence().c_str());
	fprintf(file, "%s", ";RE=");
	fprintf(file, "%i", SV->get_support());
	fprintf(file, "%s", "\tGT:DR:DV\t./.:.:");
	fprintf(file, "%i", SV->get_support());
	fprintf(file, "%c", '\n');

}

