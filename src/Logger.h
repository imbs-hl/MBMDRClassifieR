
#ifndef LOGGER_H_
#define LOGGER_H_

#include <string>
#include <sstream>
#include <mutex>
#include <memory>
#include <fstream>
#include <iostream>

// log message levels
enum Level	{ Config, Info, Warning, Severe };
class Logger;

class Logstream : public std::ostringstream {
public:
  Logstream(Logger& log, Level level);
  Logstream(const Logstream& ls);
  ~Logstream();

private:
  Logger& log;
  Level level;
};

class Logger {
public:
  Logger(std::string filename);
  Logger(std::string filename, unsigned int verboselevel);
  virtual ~Logger();

  void log(Level level, std::string message);
  void log(Level level, std::string message, unsigned int verboselevel);

  Logstream operator()();
  Logstream operator()(Level nLevel);

protected:
  std::ostream& getOutputWriter();

private:
  const tm* getLocalTime();

private:
  std::mutex mutex;
  std::ofstream stream;

  unsigned int verbose_level;

  tm local_time;
};

#endif /* LOGGER_H_ */
