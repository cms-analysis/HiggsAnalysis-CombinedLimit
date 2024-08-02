#include "../interface/CombineLogger.h"
using namespace std;

// counter for Logger calls
int CombineLogger::nLogs=0;

const char*  CombineLogger::fName = "combine_logger.out";

CombineLogger* CombineLogger::pL = nullptr;

CombineLogger& CombineLogger::instance()
{

	if (pL == nullptr)
		pL = new CombineLogger();

	return *pL;
}

CombineLogger::CombineLogger()
{
	outStream.open(fName, ios_base::out);
}

void CombineLogger::log(const std::string & _file, const int _lineN, const string& _logmsg, const string& _function)
{
	std::cout << _logmsg << std::endl;
	outStream << _file << "[" << _lineN << "] " << ": (in function: " << _function << ") - "  << _logmsg << endl;
	nLogs++; 
}

void CombineLogger::printLog()
{
	std::cout << nLogs << " log messages saved to " << fName << std::endl;
}

CombineLogger::~CombineLogger()
{
	// Combine will delete
	outStream.close();
	delete CombineLogger::pL;
	CombineLogger::pL = nullptr;
}