#include "Logger.h"

Logger::Logger() {
  loggerLevel = Info;
}

Logger::Logger(LoggerType loggerLevel_init) {
  loggerLevel = loggerLevel_init;
}

void Logger::error(std::string text) {
  if (loggerLevel >= Error) {
    std::cout << "ERROR: "<< text << std::endl;
  }
}


void Logger::info(std::string text) {
  if (loggerLevel >= Info) {
    std::cout << "INFO: " << text << std::endl;
  }
}

void Logger::debug(std::string text) {
  if (loggerLevel >= Debug) {
    std::cout << "DEBUG: " << text << std::endl;
  }
}

void Logger::verbose(std::string text) {
  if(loggerLevel >= Verbose) {
    std::cout << "VERBOSE: " << text << std::endl;
  }
}

void Logger::print(std::string text) {
  if(loggerLevel != None) {
    std::cout << text << std::endl;
  }
}
