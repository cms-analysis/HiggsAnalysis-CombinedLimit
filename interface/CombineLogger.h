#ifndef HiggsAnalysis_CombinedLimit_CombineLogger_h
#define HiggsAnalysis_CombinedLimit_CombineLogger_h

#include <iostream>
#include <fstream>
#include <string>

class CombineLogger
{
	public: 
		static int nLogs;

		static CombineLogger& instance();
		
		static void setName(const char* _fName){
			fName=_fName;
		};

		void log(const std::string & _file, const int _lineN, const std::string& _logmsg, const std::string& _function);
		void printLog();

	protected:
		// Static variable for the instance  
		static CombineLogger* pL;

		static const char*  fName;
		std::ofstream outStream;
		CombineLogger();
		virtual ~CombineLogger();
};
#endif