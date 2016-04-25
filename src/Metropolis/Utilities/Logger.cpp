#include "Logger.h"

Logger::Logger() {
	loggerLevel = NONE;
}

Logger::Logger(LoggerType loggerLevel_init) {
	loggerLevel = loggerLevel_init;
}

void Logger::verbose(std::string text) {
	if(loggerLevel == VERBOSE) {
		std::cout << text << std::endl;
	}
}

void Logger::print(std::string text) {
	if(loggerLevel != NONE) {
		std::cout << text << std::endl;
	}
}