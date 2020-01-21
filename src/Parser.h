/*
 * Parser.h
 *
 *  Created on: May 29, 2012
 *      Author: fritz
 */

#ifndef PARSER_H_
#define PARSER_H_

#include "Alignment.h"

struct str_read{
	uint chr_id;
	uint start;
	ushort length;
	ushort strand; //1=forward, 2=backward;
};

class Parser {
public:
	virtual ~Parser(){};
	virtual Alignment * parseRead(uint16_t mappingQv) = 0;
	virtual void parseReadFast(uint16_t mappingQv,Alignment *& align) = 0;
};

#endif /* PARSER_H_ */
