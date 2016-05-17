//#include <string>
//#include <sstream>

#include "Parsing.h"

std::string getExtension(const std::string& filepath)
{
	unsigned last_dot = filepath.find_last_of(".");
	unsigned last_delim = filepath.find_last_of("/\\");
	std::string extension;
	if (filepath.find_last_of("/\\") == std::string::npos || last_delim < last_dot)
		extension = filepath.substr(last_dot + 1);
	return extension;
}

template<>
bool fromString<bool>(const std::string& str, bool& result)
{
	if (str == "True" || str == "true" || str == "TRUE" || str == "1")
	{
		result = true;
		return true;
	}
	else if (str == "False" || str == "false" || str == "FALSE" || str == "0")
	{
		result = false;
		return true;
	}
	else
	{
		return false;
	}
}

template <>
bool fromString<std::string>(const std::string& str, std::string& result)
{
	result = str;
	return true;
}
