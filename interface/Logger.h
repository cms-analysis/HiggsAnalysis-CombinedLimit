/* 
 * Logger adapted from 
 * Professional C++, 2nd Edition, Oct 2011
 * Marc Gregoire, Nicholas A. Solter, Scott J. Kleper
 * ISBN: 978-0-470-93244-5
 * http://www.wiley.com/WileyCDA/WileyTitle/productCd-0470932449.html
*/

// Logger.h
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mutex>

class Logger
{
public:
	static const std::string kLogLevelDebug;
	static const std::string kLogLevelInfo;
	static const std::string kLogLevelError;

	static int nLogLevelInfo;
	static int nLogLevelDebug;
	static int nLogLevelError;

	// Returns a reference to the singleton Logger object
	static Logger& instance();

	// Logs a single message at the given log level
	void log(const std::string& inMessage, 
		const std::string& inLogLevel,
		const std::string& inFunction);

	// Logs a vector of messages at the given log level
	void log(const std::vector<std::string>& inMessages, 
		const std::string& inLogLevel,
		const std::string& inFunction);
	
	void printLog();

protected:
	// Static variable for the one-and-only instance  
	static Logger* pInstance;

	// Constant for the filename
	static const char* const kLogFileName;

	// Data member for the output stream
	std::ofstream mOutputStream;

	// Embedded class to make sure the single Logger
	// instance gets deleted on program shutdown.
	friend class Cleanup;
	class Cleanup
	{
	public:
		~Cleanup();
	};

	// Logs message. The thread should own a lock on sMutex
	// before calling this function.
	void logHelper(const std::string& inMessage, 
		const std::string& inLogLevel,
		const std::string& inFunction);
	

private:
	Logger();
	virtual ~Logger();
	Logger(const Logger&);
	Logger& operator=(const Logger&);
	static std::mutex sMutex;
};
