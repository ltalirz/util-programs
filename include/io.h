/*
 * io.h
 */
#ifndef IO_H
#define IO_H

#include <string>

namespace io {

typedef std::string String;
typedef std::vector<char> Binary;

bool readFile(String filename, String& content);

bool readBinary(String filename, Binary& content);



}
#endif
