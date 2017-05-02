#include <iostream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <Rcpp.h>
#include "Logger.h"

// needed for MSVC
#ifdef WIN32
#define localtime_r(_Time, _Tm) localtime_s(_Tm, _Time)
#endif // localtime_r

// Convert date and time info from tm to a character string
// in format "YYYY-mm-DD HH:MM:SS" and send it to a stream
std::ostream& operator<< (std::ostream& stream, const tm* tm) {
	// I had to muck around this section since GCC 4.8.1 did not implement std::put_time
	//	return stream << std::put_time(tm, "%Y-%m-%d %H:%M:%S");
	return stream << 1900 + tm->tm_year << '-' <<
			std::setfill('0') << std::setw(2) << tm->tm_mon + 1 << '-'
			<< std::setfill('0') << std::setw(2) << tm->tm_mday << ' '
			<< std::setfill('0') << std::setw(2) << tm->tm_hour << ':'
			<< std::setfill('0') << std::setw(2) << tm->tm_min << ':'
			<< std::setfill('0') << std::setw(2) << tm->tm_sec;
}

Logger::Logger(std::string filename) : verbose_level(1) {
	if(filename.size()) {
		stream.open(filename, std::fstream::out | std::fstream::app | std::fstream::ate);
	}
}

Logger::Logger(std::string filename, unsigned int verbose_level) : verbose_level(verbose_level) {
	if(filename.size()) {
		stream.open(filename, std::fstream::out | std::fstream::app | std::fstream::ate);
	}
}

Logger::~Logger() {
	if(stream) {
		stream.flush();
		stream.close();
	}
}

Logstream Logger::operator()() {
	return Logstream(*this, Info);
}

Logstream Logger::operator()(Level nLevel) {
	return Logstream(*this, nLevel);
}

const tm* Logger::getLocalTime() {
	auto in_time_t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	localtime_r(&in_time_t, &local_time);
	return &local_time;
}

void Logger::log(Level level, std::string message) {
	const static char* LevelStr[] = { "Config", "Info", "Warning", "Severe" };

	mutex.lock();
	getOutputWriter() << '[' << getLocalTime() << ']'
			<< '[' << LevelStr[level] << "]\t"
			<< message << std::endl;
	mutex.unlock();
}

void Logger::log(Level level, std::string message, unsigned int verbose_level) {
	const static char* LevelStr[] = { "Config", "Info", "Warning", "Severe" };

	if(this->verbose_level >= verbose_level) {
		mutex.lock();
		getOutputWriter() << '[' << getLocalTime() << ']'
				<< '[' << LevelStr[level] << "]\t"
				<< message << std::endl;
		mutex.unlock();
	}
}

std::ostream& Logger::getOutputWriter() {
	if (stream.is_open()) {
		return stream;
	} else {
		return std::cout;
	}
}

Logstream::Logstream(Logger& log, Level level) :
				log(log), level(level) { }

Logstream::Logstream(const Logstream& ls) :
				log(ls.log), level(ls.level) {
	// As of GCC 8.4.1 basic_stream is still lacking a copy constructor
	// (part of C++11 specification)
	//
	// GCC compiler expects the copy constructor even thought because of
	// RVO this constructor is never used
}

Logstream::~Logstream() {
	log.log(level, str());
}
