#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <boost/lexical_cast.hpp>

using namespace std;

	
	struct dat {
	  string key;
	  string value;
	};

	class FileReader
	{
		public:
			FileReader()
			{

			}
			void readParameters(string filename)
			{
				ifstream input;
				string line;

				input.open(filename, ios::in);
				if (input.is_open())
				{
					while (getline(input, line))
					{
						size_t space = line.find_first_of(' ');
						dat param;
						param.key = line.substr(0, space);
						param.value = line.substr(space + 1);
						m_params.push_back(param);
						//cout << param.key << " " << param.value << '\n';
					}
					input.close();
				}
			}
			
			template<typename Type>
			Type getParameter(string valueName)
			{
				for(vector<dat>::iterator it = m_params.begin() ; it != m_params.end(); ++it)
				{
					if(!it->key.compare(valueName))
					{
						return boost::lexical_cast<Type>(it->value);
					}
				}
				return boost::lexical_cast<Type>("nix wars");
			}
		
		private:
			vector<dat> m_params;
	};


