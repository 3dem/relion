/***************************************************************************
 *
 * Author: "Dari Kimanius"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/


#ifndef PYTHON_DEPENDENCIES_H
#define PYTHON_DEPENDENCIES_H

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <iostream>

// Helper macro for stringification
#define STRINGIFY(x) #x

namespace python_dependencies
{
	std::string get_python_exe_path()
	{
		// First check if the environmental variable is set
		const char* python_exe_path = getenv("RELION_PYTHON_EXECUTABLE");
		if (python_exe_path != nullptr)
			return std::string{python_exe_path};

		// Then check if pre-compiled variables are set
#ifdef PYTHON_EXE_PATH
#define PYTHON_EXE_PATH_STRING STRINGIFY(PYTHON_EXE_PATH)
		return PYTHON_EXE_PATH_STRING;
#else
		// If not, return default
		return "python";
#endif
	}

	std::string get_torch_home_path()
	{
		// First check if the environmental variable is set
		const char* torch_home_path = getenv("TORCH_HOME");
		if (torch_home_path != nullptr)
			return std::string{torch_home_path};

		// Then check if pre-compiled variables are set
#ifdef TORCH_HOME_PATH
#define TORCH_HOME_PATH_STRING STRINGIFY(PYTHON_EXE_PATH)
		return TORCH_HOME_PATH_STRING;
#else
		// If not, return empty string
		return "";
#endif
	}

	class SetEnvException : public std::exception
	{
	private:
		std::string variableName_;
		std::string variableValue_;
		mutable std::string errorMessage_;

	public:
		SetEnvException(std::string  variableName, std::string  variableValue)
				: variableName_(std::move(variableName)), variableValue_(std::move(variableValue)) {}

		const char* what() const noexcept override {
			errorMessage_ = "Failed to set environmental variable: " + variableName_ + " = " + variableValue_;
			return errorMessage_.c_str();
		}
	};

	class InterpreterException : public std::exception
	{
	private:
		std::string cmd_;
		mutable std::string errorMessage_;

	public:
		explicit InterpreterException(std::string cmd="")
				: cmd_(std::move(cmd)) {}

		const char* what() const noexcept override {
			errorMessage_ =
			"---------------------------------- PYTHON ERROR ---------------------------------\n"
			"   Has RELION been provided a Python interpreter with the correct environment?   \n"
			" The interpreter can be passed to RELION either during Cmake configuration with  \n"
			"     using the Cmake flag -DPYTHON_EXE_PATH=<path to python interpreter> or      \n"
			"          by setting the environmental variable RELION_PYTHON_EXECUTABLE.        \n"
			"---------------------------------------------------------------------------------\n"
			"Failed to execute command: " + cmd_;
			return errorMessage_.c_str();
		}
	};

	void export_torch_home_path()
	{
		std::basic_string<char> torch_home_path = get_torch_home_path();
		if (setenv("TORCH_HOME", torch_home_path.c_str(), 1) != 0)
			throw SetEnvException("TORCH_HOME", torch_home_path);
	}

	std::string get_full_cmd(const std::string& python_cmd)
	{
		std::string python_exe_path = get_python_exe_path();
		std::string torch_home_path = get_torch_home_path();
		return "TORCH_HOME=" + torch_home_path + " " + python_exe_path + " " + python_cmd;
	}

	std::string execute(const std::string& python_cmd)
	{
		std::string cmd = get_full_cmd(python_cmd);

		FILE* pipe = popen(cmd.c_str(), "r");
		if (!pipe)
			throw std::runtime_error("Failed to dispatch command: " + cmd);

		char buffer[128];
		std::string result = "";

		// read till end of process:
		while (!feof(pipe))
			if (fgets(buffer, 128, pipe) != NULL)
				result += buffer;

		pclose(pipe);

		return result;
	}

}
#endif //PYTHON_DEPENDENCIES_H