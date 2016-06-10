#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <string>

enum LoggerType {
  //VERBOSE, NON_VERBOSE, NONE  
  None,
  Error,
  Info,
  Debug,
  Verbose
};  

class Logger {
  private:
    /**
     * This stores the level at which the logger prints.
     */
    LoggerType loggerLevel;

  public:

    /**
     * Default constructor for logger.
     */
    Logger();

    /**
     * Constructor for the Logger class. Initializes the loggerLevel.
     * @param loggerLevel_init The loggerType to initialize the logger to.
     */
    Logger(LoggerType loggerLevel_init);

    /**
     * Logs error messages
     */
    void error(std::string text);

    /**
     * Logs informational messages
     */
    void info(std::string text);

    /**
     * Logs debug messages
     */
    void debug(std::string text);

    /**
     * Prints messages that only a verbose level logger will print.
     */
    void verbose(std::string text);

    /**
     * Prints messages that only a verbose level or non verbose level logger
     * will print.
     */
    void print(std::string text);
};


#endif
