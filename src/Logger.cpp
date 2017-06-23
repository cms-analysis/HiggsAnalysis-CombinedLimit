/* 
 * Logger adapted from 
 * Professional C++, 2nd Edition, Oct 2011
 * Marc Gregoire, Nicholas A. Solter, Scott J. Kleper
 * ISBN: 978-0-470-93244-5
 * http://www.wiley.com/WileyCDA/WileyTitle/productCd-0470932449.html
*/

#include <stdexcept>
#include "HiggsAnalysis/CombinedLimit/interface/Logger.h"

using namespace std;

const string Logger::kLogLevelDebug = "DEBUG";
const string Logger::kLogLevelInfo = "INFO";
const string Logger::kLogLevelError = "ERROR";

int Logger::nLogLevelInfo=0;
int Logger::nLogLevelDebug=0;
int Logger::nLogLevelError=0;

const char* const Logger::kLogFileName = "combine_logger.out";

Logger* Logger::pInstance = nullptr;

mutex Logger::sMutex;

Logger& Logger::instance()
{
	static Cleanup cleanup;

	lock_guard<mutex> guard(sMutex);
	if (pInstance == nullptr)
		pInstance = new Logger();

	return *pInstance;
}

Logger::Cleanup::~Cleanup()
{
	lock_guard<mutex> guard(Logger::sMutex);
	delete Logger::pInstance;
	Logger::pInstance = nullptr;
}

Logger::~Logger()
{
	mOutputStream.close();
}

Logger::Logger()
{
	mOutputStream.open(kLogFileName, ios_base::out);
	if (!mOutputStream.good()) {
		throw runtime_error("Unable to initialize the Logger!");
	} 
}

void Logger::log(const string& inMessage, const string& inLogLevel, const string& inFunction)
{
	lock_guard<mutex> guard(sMutex);
	logHelper(inMessage, inLogLevel, inFunction);
}

void Logger::log(const vector<string>& inMessages, const string& inLogLevel, const string& inFunction)
{
	lock_guard<mutex> guard(sMutex);
	for (size_t i = 0; i < inMessages.size(); i++) {
		logHelper(inMessages[i], inLogLevel, inFunction);
	}
}

void Logger::logHelper(const std::string& inMessage, const std::string& inLogLevel, const std::string& inFunction)
{
	mOutputStream << inLogLevel << ": (function: " << inFunction << ") "  << inMessage << endl;
	if (inLogLevel == kLogLevelInfo)  nLogLevelInfo++;
	if (inLogLevel == kLogLevelDebug) nLogLevelDebug++;
	if (inLogLevel == kLogLevelError) nLogLevelError++;
}

void Logger::printLog()
{
	std::cout << "Printing Message Summary From ... " << kLogFileName << std::endl;
	std::cout << "----------------------------------------------" << std::endl;
	std::cout << "Messages of type " << kLogLevelInfo  << " : " << nLogLevelInfo  << std::endl;
	std::cout << "Messages of type " << kLogLevelDebug << " : " << nLogLevelDebug << std::endl;
	std::cout << "Messages of type " << kLogLevelError << " : " << nLogLevelError << std::endl;
	std::cout << "----------------------------------------------" << std::endl;
	std::ifstream f(kLogFileName);
	if (f.is_open()) std::cout << f.rdbuf();
}
