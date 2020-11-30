#ifndef DYN_LOG_H
#define DYN_LOG_H

#include <string>
#include <sys/time.h>



class Logger;
class IOParser;

class Log
{
	public:
		
		static void print(std::string s);
		static void extend(std::string s);
		static void warn(std::string s);
		
		static void printBinaryChoice(
				std::string start, bool condition, 
				std::string endIfTrue, std::string endIfFalse);
		
		static void beginSection(std::string name);
		static void endSection();
		
		static void beginProgress(std::string name, long int totalItems);
		static void updateProgress(long int itemsDone);
		static void endProgress();
		
		static void readParams(IOParser& parser);
		static void saveSettingsFile(
				const std::string& contents,
				const std::string& directory,
				const std::string& file);
		
		
	private:
		
		static Logger* logger;
		static long int lastItemsDone_global, totalItems_global;
};

class Logger
{
	public:
		
		virtual void print(std::string s) = 0;
		virtual void extend(std::string s) = 0;
		virtual void warn(std::string s) = 0;
		
		virtual void beginSection(std::string name) = 0;
		virtual void endSection() = 0;
		
		virtual void beginProgress(std::string name, long int totalItems) = 0;
		virtual void updateProgress(long int itemsDone) = 0;
		virtual void endProgress() = 0;
};

class MinimalistLogger : public Logger
{
	public:
		
		MinimalistLogger();
		
		
		void print(std::string s);
		void extend(std::string s);
		void warn(std::string s);
		
		void beginSection(std::string name);
		void endSection();
		
		void beginProgress(std::string name, long int totalItems);
		void updateProgress(long int itemsDone);
		void endProgress();
		
		
	private:
		
		int level, tabSize;
		long int totalItems;
		std::string head;
};

class FancyLogger : public Logger
{
	public:
		
		FancyLogger();
		
		
		void print(std::string s);
		void extend(std::string s);
		void warn(std::string s);
		
		void beginSection(std::string name);
		void endSection();
		
		void beginProgress(std::string name, long int totalItems);
		void updateProgress(long int itemsDone);
		void endProgress();
		
		
		void setColours(std::string tree, std::string accent);
		
		
	private:
		
		long int totalItems;
		int barWidth, indentLevel;
		std::string head, textColour, treeColour, accentColour, dotStr;
		
		time_t startt, prevt;
		
		
		const std::string 
			dot = "■", block = "■";
		
		const std::string 
		
			black   = "\u001b[30m",
			red     = "\u001b[31m",
			green   = "\u001b[32m",
			yellow  = "\u001b[33m",
			blue    = "\u001b[34m",
			magenta = "\u001b[35m",
			cyan    = "\u001b[36m",
			white   = "\u001b[37m",
		
			reset = "\u001b[0m",
		
			bright_black   = "\033[90m",
			bright_red     = "\033[91m",
			bright_green   = "\033[92m",
			bright_yellow  = "\033[93m",
			bright_blue    = "\033[94m",
			bright_magenta = "\033[95m",
			bright_cyan    = "\033[96m",
			bright_white   = "\033[97m";
		
		
		std::string entryChar();
		std::string continueChar();
		
};

#endif
