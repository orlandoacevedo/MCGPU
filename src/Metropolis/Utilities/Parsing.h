#ifndef METROSIM_PARSING_H
#define METROSIM_PARSING_H

#include <string>
#include <sstream>

std::string getExtension(const std::string& filepath);

template <typename T>
bool fromString(const std::string& str, T& result)
{
	std::istringstream stream(str);
	if (!(stream >> result))
		return false;
	return true;
}

template <typename T>
bool fromString(char* cstr, T& result)
{
	if (!cstr)
		return false;
	std::string str(cstr);
	return fromString<T>(str, result);
}

template<>
bool fromString<bool>(const std::string& str, bool& result);

template <>
bool fromString<std::string>(const std::string& str, std::string& result);

template <typename T>
bool toString(T& val, std::string& result)
{
	std::ostringstream stream;
	if (!(stream << val))
		return false;
	result = stream.str();
	return true;
}

#endif
