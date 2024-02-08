/*
this header contain function prototypes that helps to parse a .ini file using STL map container and template functions 
*/
#ifndef INIPARSER_H
#define INIPARSER_H

#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

class INIParser 
{
    public:
    static std::unordered_map<std::string, std::unordered_map<std::string, std::string>> parse(const std::string& filename);
    static int getInt(const std::unordered_map<std::string, std::string>& section, const std::string& key);
    static double getDouble(const std::unordered_map<std::string, std::string>& section, const std::string& key);
    static std::string getString(const std::unordered_map<std::string, std::string>& section, const std::string& key);

    private:
    static void trim(std::string& str);
};

#endif // INIPARSER_H
