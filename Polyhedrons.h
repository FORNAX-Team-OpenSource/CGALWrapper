#pragma once

#include "Polyhedron.h"
#include <list>
#include <vector>
#include <string>

namespace JSONModels
{
	class Polyhedrons : public JSONBase
	{
	public:
		virtual ~Polyhedrons() {};
		virtual bool Deserialize(const std::string& s);

		// Getters/Setters.
		std::vector<Polyhedron>& PolyhedronList() { return _polyhedrons; }
	public:
		virtual bool Deserialize(const rapidjson::Value& obj) { return false; };
		virtual bool Serialize(rapidjson::Writer<rapidjson::StringBuffer>* writer) const;
	private:
		std::vector<Polyhedron> _polyhedrons;
	};
}

