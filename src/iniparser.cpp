#include "iniparser.h"


std::unordered_map<std::string, std::unordered_map<std::string, std::string>> INIParser::parse(const std::string& filename) 
{
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> sections;
    std::string currentSection;

    std::string line;
    while (std::getline(file, line)) 
    {
        if (line.empty() || line[0] == ';' || line[0] == '#') 
        {
            continue;
        }

        if (line[0] == '[' && line[line.length() - 1] == ']') 
        {
            currentSection = line.substr(1, line.length() - 2);
            continue;
        }

        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, '=') && std::getline(iss, value)) 
        {
            trim(key);
            trim(value);
            sections[currentSection][key] = value;
        }
    }

    return sections;
}

int INIParser::getInt(const std::unordered_map<std::string, std::string>& section, const std::string& key) 
{
    return std::stoi(section.at(key));
}

double INIParser::getDouble(const std::unordered_map<std::string, std::string>& section, const std::string& key) 
{
    return std::stod(section.at(key));
}

std::string INIParser::getString(const std::unordered_map<std::string, std::string>& section, const std::string& key) 
{
    return section.at(key);
}

void INIParser::trim(std::string& str) 
{
    size_t first = str.find_first_not_of(" \t");
    size_t last = str.find_last_not_of(" \t");
    str = str.substr(first, (last - first + 1));
}
