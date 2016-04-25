#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <string>

enum LoggerType {
	VERBOSE, NON_VERBOSE, NONE	
};	

class Logger {
	private:
		/**
		 * This stores the level at which the logger prints.
		 */
		LoggerType loggerLevel;
		
	public:
		
		/**
		 * Blank constructor for logger.
		 */
		Logger();
		
		/**
		 * Constructor for the Logger class. Initializes the loggerLevel.
		 * @param loggerLevel_init The loggerType to initialize the logger to.
		 */
		Logger(LoggerType loggerLevel_init);
		
		/**
		 * Prints messages that only a verbose level logger will print.
		 */
		void verbose(std::string text);
		
		/**
		 * Prints messages that only a verbose level or non verbose level logger will print.
		 */
		void print(std::string text);
};


#endif
