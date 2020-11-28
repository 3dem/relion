#include "log.h"
#include "ansi_codes.h"
#include <src/time.h>
#include <src/args.h>
#include <iomanip>

Logger* Log::logger = new FancyLogger();
long int Log::lastItemsDone_global = -1;
long int Log::totalItems_global = 0;

#define MAX_UPDATES_PER_BAR 200


void Log::print(std::string s)
{
	logger->print(s);
}

void Log::extend(std::string s)
{
	logger->extend(s);
}

void Log::warn(std::string s)
{
	logger->warn(s);
}

void Log::printBinaryChoice(std::string start, bool condition, std::string endIfTrue, std::string endIfFalse)
{
	if (condition)
	{
		logger->print(start + endIfTrue);
	}
	else
	{
		logger->print(start + endIfFalse);
	}
}

void Log::beginSection(std::string name)
{
	logger->beginSection(name);
}

void Log::endSection()
{
	logger->endSection();
}

void Log::beginProgress(std::string name, long int totalItems)
{
	lastItemsDone_global = -1;
	totalItems_global = totalItems;
	
	logger->beginProgress(name, totalItems);
}

void Log::updateProgress(long int itemsDone)
{
	if ( lastItemsDone_global > 0 && 
	    ((itemsDone - lastItemsDone_global) < (totalItems_global / MAX_UPDATES_PER_BAR))) 
	{
		return;
	}
	
	lastItemsDone_global = itemsDone;
	
	logger->updateProgress(itemsDone);
}

void Log::endProgress()
{
	lastItemsDone_global = -1;
	totalItems_global = 0;
	
	logger->endProgress();
}

void Log::readParams(IOParser &parser)
{
	int gen_section = parser.addSection("Display options");
	
	{
		using namespace std;
		
		vector<pair<string, pair<string, string>>> choices(0);
		choices.reserve(16);
		
		choices.push_back(make_pair("molten", make_pair(ANSI_BLACK, ANSI_ORANGE_2_256)));
		choices.push_back(make_pair("metal", make_pair(ANSI_GREY_5_256, ANSI_GREY_11_256)));
		choices.push_back(make_pair("starlight", make_pair(ANSI_CLEAN_BLUE_1_256, ANSI_SKY_BLUE_1_256)));
		choices.push_back(make_pair("inferno", make_pair(ANSI_RED, ANSI_ORANGE_1_256)));
		choices.push_back(make_pair("gfp", make_pair(ANSI_GREEN, ANSI_BRIGHT_GREEN_256)));
		choices.push_back(make_pair("mcherry", make_pair(ANSI_DEEPEST_PURPLE_256, ANSI_PURPLE_256)));
		choices.push_back(make_pair("hal", make_pair(ANSI_BLACK, ANSI_FULL_RED_256)));
		choices.push_back(make_pair("sal", make_pair(ANSI_BLACK, ANSI_CLEAN_BLUE_4_256)));
		choices.push_back(make_pair("summer", make_pair(ANSI_SKY_BLUE_1_256, ANSI_ORANGE_4_256)));
		choices.push_back(make_pair("winter", make_pair(ANSI_SKY_BLUE_3_256, ANSI_PALE_YELLOW_2_256)));
		choices.push_back(make_pair("monochrome", make_pair(ANSI_RESET, ANSI_RESET)));
		
		std::string homeDir = std::getenv("HOME");
		std::string settingsDir = homeDir+"/.relion";
		std::string settingsFile = settingsDir + "/theme.conf";
		
		std::string defaultTheme = choices[0].first;
		
		std::ifstream file(settingsFile);
		if (file)
		{
			file >> defaultTheme;
		}
	
		
		std::string themeStr = parser.getOption("--theme", "Select display theme", defaultTheme);
		bool do_list_themes = parser.checkOption("--list_themes", "List available themes and quit");
		
		if (do_list_themes)
		{
			std::cout << "\nThe following themes are available:\n";
			
			FancyLogger fl;
			
			for (int i = 0; i < choices.size(); i++)
			{
				fl.setColours(choices[i].second.first, choices[i].second.second);
				
				fl.beginSection(choices[i].first);
				fl.endSection();
			}
			
			std::cout << '\n';
			
			MinimalistLogger ml;
			
			ml.beginSection("classic");
			ml.endSection();
			
			std::cout << "\n\n"
			<< "Note: \n\n"
			<< "Your choice will be stored in " << settingsFile << ".\n"
			<< "The colours will look different on different terminals.\n" << std::endl;
			
			std::exit(0);
		}
		
		bool understood = false;
		
		if (themeStr == "classic")
		{
			delete logger;
			logger = new MinimalistLogger;
			
			understood = true;
		}
		
		for (int i = 0; i < choices.size(); i++)
		{
			if (themeStr == choices[i].first)
			{
				FancyLogger* fl = new FancyLogger;
				fl->setColours(choices[i].second.first, choices[i].second.second);
						
				delete logger;
				logger = fl;
				
				understood = true;
				break;
			}
		}
		
		if (!understood)
		{
			FancyLogger* fl = new FancyLogger;
			
			delete logger;
			logger = fl;
			
			fl->setColours(choices[0].second.first, choices[0].second.second);
			
			
			if (themeStr == defaultTheme)
			{
				logger->warn("Unknown display theme " + themeStr + " found in " + settingsFile);
				logger->extend("Reverting to '" + choices[0].first + "'");

				saveSettingsFile(choices[0].first, settingsDir, settingsFile);
			}
			else
			{
				logger->print("Unknown display theme: " + themeStr);
				logger->extend("Reverting to '" + choices[0].first + "'");
			}
		}
		else
		{
			if (defaultTheme != themeStr)
			{
				logger->extend("");
				logger->print("Saving theme choice '" + std::string(ANSI_BOLD) + themeStr 
							  + std::string(ANSI_RESET) + "' to " + settingsFile);
				
				saveSettingsFile(themeStr, settingsDir, settingsFile);
			}
		}
	}
	
	std::cout << '\n';
}

void Log::saveSettingsFile(const std::string &contents, const std::string &directory, const std::string &file)
{
	int res0 = system((std::string("mkdir -p ")+directory).c_str());

	if (res0)
	{
		logger->warn("Unable to create UI settings directory " + directory);
	}
	else
	{
		int res1 = system((std::string("echo ")+contents+" > "+file).c_str());

		if (res1)
		{
			logger->warn("Unable to write to UI settings file " + file);
		}
	}
}



MinimalistLogger::MinimalistLogger()
:	level(0), tabSize(3), totalItems(0),
	head(" + ")
{
}

void MinimalistLogger::print(std::string s)
{
	std::cout << head << s << std::endl;
}

void MinimalistLogger::extend(std::string s)
{
	std::cout << std::string(head.length(), ' ') << s << std::endl;
}

void MinimalistLogger::warn(std::string s)
{
	std::cout << std::string(head.length(), ' ') << '!' << s << std::endl;
}

void MinimalistLogger::beginSection(std::string name)
{
	std::cout << head << name << ':' << std::endl;
	
	level++;
	
	head = std::string(level * tabSize, ' ');
	
	if (level == 0) head = head + " + ";
	else if (level == 1) head = head + " - ";
	else head = head + "   ";
}

void MinimalistLogger::endSection()
{
	level--;
	
	head = std::string(level * tabSize, ' ');
	
	if (level == 0) head = head + " + ";
	else if (level == 1) head = head + " - ";
	else head = head + "   ";
}

void MinimalistLogger::beginProgress(std::string name, long int totalItems)
{
	std::cout << head << name << ':' << std::endl;
	init_progress_bar(totalItems);
	
	this->totalItems = totalItems;
}

void MinimalistLogger::updateProgress(long int itemsDone)
{
	progress_bar(itemsDone);
}

void MinimalistLogger::endProgress()
{
	progress_bar(totalItems);
}


/*
  
  Prototype:
  
  
■ initialising                 // beginSection(0, "initialising")
├──────────────┘
│
├■ component X                 // beginSection(1, "component X")
│├─────────────┘
│├■ loading X                  // print("loading X")
│├■ transmogrifying            // beginProgress("transmogrifying", N)

│││ [         ] 0%,                     
│││ [         ] 0%,                     // after updateProgress(0)
│││ [■■■■■■   ] 56%,  41 sec remaining  // after updateProgress(i)
│││ [■■■■■■■■■] 100%                    // after updateProgress(N-1) or endProgress()

││└                           // endProgress()  
│└                            // endSection(1)
├■ component Y                 // beginSection(1, "component Y")
│├─────────────┘
│├■ loading Y
││
│├■ analysing Y                // beginSection(2, "analysing Y")
││├─────────────┘
││├■ test0
││├■ test1
││└            
│└            
└


Alternatives:


■ initialising                 // beginSection(0, "initialising")
├═─═─═─═─═─═─═─═─
├■ component X                 // beginSection(1, "component X")
│├═─═─═─═─═─═─═─═─
│├■ loading X                  // print("loading X")
│├■ transmogrifying            // beginProgress("transmogrifying", N)


■ initialising                 // beginSection(0, "initialising")
╠────────────────
╠■ component X                 // beginSection(1, "component X")
║╠───────────────
║╠■ loading X                  // print("loading X")
║╠■ transmogrifying            // beginProgress("transmogrifying", N)
║



 label:
[■■■■■■   ]

or 

 label:───┐
└██████░░░┘

or

 label: 
║██████░░░║

 label:
└▀▀▀▀▀    ┘


■ initialising         
╔══════════════╗
├■ component X ║    
│╔══════════════╗
│├■ part 1      ║
│└──────────────╝
│╔══════════════╗
│├■ part 2      ║
│└──────────────╝
└──────────────┘


■ initialising         
┌──────────────╗
├■ component X    
│┌─────────╗
│├■ part 1 │     
│╚─────────┘
│┌─────────╗
│├■ part 2 │     
│╚─────────┘
╚──────────────┘




│ ┤ ┐ └ ┴ ┬ ├ ─ ┼ ┘┌

╣ ║ ╗ ╝ ╚ ╔ ╩ ╦ ╠ ═ ╬

╡ ╢ ╖ ╕ ╜ ╛ ╞ ╟ ╧ ╨ ╤ ╥ ╙ ╘ ╒ ╓ ╫ ╪


░ ▒ ▓ │ ┤ ╡ ╢ ╖ ╕ ╣ ║ ╗ ╝ ╜ ╛ ┐

└ ┴ ┬ ├ ─ ┼ ╞ ╟ ╚ ╔ ╩ ╦ ╠ ═ ╬ ╧

╨ ╤ ╥ ╙ ╘ ╒ ╓ ╫ ╪ ┘ ┌ █ ▄ ▌ ▐ ▀

*/

FancyLogger::FancyLogger()
:	totalItems(0),
	barWidth(40),
	indentLevel(0),
	head("")
{
	textColour = reset;		
	treeColour = ANSI_BLACK;
	accentColour = ANSI_BLUE;
	
	head = treeColour;	
	dotStr = accentColour + dot;
}

void FancyLogger::print(std::string s)
{
	std::cout << head << entryChar() << dotStr << " " << textColour << s << std::endl;
}

void FancyLogger::extend(std::string s)
{
	std::cout << head << continueChar() << "  " << textColour << s << std::endl;
}

void FancyLogger::warn(std::string s)
{
	std::cout << head << continueChar() << accentColour << "! " << textColour << s << std::endl;
}

void FancyLogger::beginSection(std::string name)
{
	std::cout 
		<< head << continueChar() << '\n' 
		<< head << entryChar() << dotStr << " " << textColour << name << "\n";
	
	indentLevel++;
	
	head = treeColour + " ";
	
	for (int i = 0; i < indentLevel-1; i++)
	{
		head = head + "│";
	}
	
	std::cout << head << "╞═";
	
	for (int i = 0; i < name.length(); i++)
	{
		std::cout << "═";
	}
	
	std::cout << "═" << ANSI_RESET << std::endl;
}

void FancyLogger::endSection()
{
	indentLevel--;
	
	std::cout << head << "└\n" << ANSI_RESET;
	
	head = treeColour;
	
	if (indentLevel > 0) head = head + " ";
	
	for (int i = 0; i < indentLevel-1; i++)
	{
		head = head + "│";
	}
}

void FancyLogger::beginProgress(std::string name, long int totalItems)
{
	print(name+":");
	
	time_t currt = time(NULL);
	
	prevt = currt;
	startt = currt;
	
	this->totalItems = totalItems;
	
	std::cout << head << continueChar() << ANSI_BOLD << " [" << ANSI_RESET;
	
	for (int i = 0; i < barWidth; i++)
	{
		std::cout << ' ';
	}
		
	std::cout << treeColour << ANSI_BOLD << "] " << ANSI_RESET << textColour << "remaining time unknown.";
	std::cout << ANSI_RESET << '\r';
	std::cout.flush();
}

void FancyLogger::updateProgress(long int itemsDone)
{
	std::cout << std::string(80, ' ');
	std::cout.flush();
	
	std::cout << '\r';
	std::cout << head << continueChar() << ANSI_BOLD << " [" << ANSI_RESET << accentColour;
	
	if (totalItems < 1)
	{
		std::cout << " BAR NOT INITIALISED ";
		return;
	}
	
	int fillTo = (int)(barWidth * itemsDone / (double)(totalItems) + 0.5);
	
	if (fillTo > barWidth) fillTo = barWidth;
	
	for (int i = 0; i < fillTo; i++)
	{
		std::cout << block;
	}
	
	for (int i = fillTo; i < barWidth; i++)
	{
		std::cout << ' ';
	}
	
	std::cout << treeColour << ANSI_BOLD << "] " << ANSI_RESET << textColour;
	
	if (itemsDone > 0)
	{
		time_t currt = time(NULL);
		
		long t1 = currt - startt; // Elapsed time
		long t2 = (long)(t1 * (float)totalItems / itemsDone); // Total time
		
		int hour = 0;
		int min = 0;
		
		float h1, h2, m1, m2;
		
		if (t2 > 60)
		{
			m1 = (float)t1 / 60.0;
			m2 = (float)t2 / 60.0;
			min = 1;
			
			if (m2 > 60)
			{
				h1 = (float)m1 / 60.0;
				h2 = (float)m2 / 60.0;
				hour = 1;
				min = 0;
			}
			else hour = 0;
		}
		else min = 0;
	
					 
		if (hour)
		{
			std::cout << std::fixed << std::setfill(' ') << std::setw(3) << std::setprecision(2) << h1;
			std::cout << '/';
			std::cout << std::fixed << std::setfill(' ') << std::setw(3) << std::setprecision(2) << h2;
			std::cout << " hrs ";
		}	
		else if (min)
		{
			
			std::cout << std::fixed << std::setfill(' ') << std::setw(3) << std::setprecision(2) << m1;
			std::cout << '/';
			std::cout << std::fixed << std::setfill(' ') << std::setw(3) << std::setprecision(2) << m2;
			std::cout << " min ";
		}
		else
		{
			
			std::cout << std::setfill(' ') << std::setw(4) << t1;
			std::cout << '/';
			std::cout << std::setfill(' ') << std::setw(4) << t2;
			std::cout << " sec ";
		}
	}
	
	std::cout  << ANSI_RESET << '\r';
	std::cout.flush();
	
}

void FancyLogger::endProgress()
{
	std::cout << std::string(80, ' ');
	std::cout.flush();
	
	std::cout << '\r';
	std::cout << head << continueChar() << ANSI_BOLD << " [" << ANSI_RESET << treeColour;
	
	for (int i = 0; i < barWidth; i++)
	{
		std::cout << block;
	}
	
	std::cout << treeColour << ANSI_BOLD << "] " << ANSI_RESET;
	
	time_t currt = time(NULL);	
	long t = currt - startt;
	
	if (t > 3600)
	{
		std::cout << std::fixed << std::setfill(' ') << std::setw(3) << std::setprecision(2) << t/3600.0;
		std::cout << " hrs ";
	}
	else if (t > 60)
	{
		std::cout << std::fixed << std::setfill(' ') << std::setw(3) << std::setprecision(2) << t/60.0;
		std::cout << " min ";
	}
	else
	{
		std::cout << std::setfill(' ') << std::setw(3) << t;
		std::cout << " sec ";
	}
	
	
	std::cout << '\n' << head << continueChar() << ANSI_RESET << std::endl;
	
	totalItems = 0;
}

void FancyLogger::setColours(std::string tree, std::string accent)
{
	treeColour = tree;
	accentColour = accent;
	
	head = treeColour;	
	dotStr = accentColour + dot;
}

std::string FancyLogger::entryChar()
{
	if (indentLevel == 0)
	{
		return " ";
	}
	else
	{
		return "╞";
	}
}

std::string FancyLogger::continueChar()
{
	if (indentLevel == 0)
	{
		return " ";
	}
	else
	{
		return "│";
	}
}


/*
  
  More alternatives:
  
│ ┤ ┐ └ ┴ ┬ ├ ─ ┼ ┘┌

╣ ║ ╗ ╝ ╚ ╔ ╩ ╦ ╠ ═ ╬

 ■ processing
 ├─────────────
 │
 ├■ component X
 │├──────────────
 ││
 │├■ preparing
 ││├────────────
 ││├■ step 1
 ││├■ step 2
 ││└
 │└
 └

 ■ processing
 ├═════════════
 │
 ├■ component X
 │├══════════════
 ││
 │├■ preparing
 ││├════════════
 ││├■ step 1
 ││├■ step 2
 ││└
 │└
 └

 ■ processing
 ╠─────────────
 ║
 ╠■ component X
 ║╠──────────────
 ║║
 ║╠■ preparing
 ║║╠────────────
 ║║╠■ step 1
 ║║╠■ step 2
 ║║╚
 ║╚
 ╚

 
 ┌───────────────┐
 ■ initialising
 ├───────────────┘
 │┌──────────────┐
 ├■ component X
 │├──────────────┘
 │├■ loading X
 │├■ transmogrifying:
 ││ [■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■]   4 sec                          
 ││
 │├■ component X ready
 │└
 └
 ┌             ┐
 ■ processing
 ├─────────────┘
 │┌             ┐
 ├■ component X
 │├─────────────┘
 ││┌
 │├■ preparing
 ││├────────────
 ││├■ step 1
 ││├■ step 2
 ││└
 │└
 └
 
 
 
 
 ┌
 ■ initialising
 ├───────────────
 │┌
 ├■ component X
 │├──────────────
 │├■ loading X
 │├■ transmogrifying:
 ││ [■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■]   4 sec                          
 ││
 │├■ component X ready
 │└
 └
               ┐ 
 ■ processing
 ├─────────────
 │               ┐
 ├■ component X
 │├──────────────
 ││             ┐ 
 │├■ preparing
 ││├────────────
 ││├■ step 1
 ││├■ step 2
 ││└
 │└
 └

 
  */
