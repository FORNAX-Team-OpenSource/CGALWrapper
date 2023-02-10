#include "Polyhedrons.h"

namespace JSONModels
{
	bool Polyhedrons::Deserialize(const std::string& s)
	{
		rapidjson::Document doc;
		if (!InitDocument(s, doc))
			return false;

		if (!doc.IsArray())
			return false;

		for (rapidjson::Value::ConstValueIterator itr = doc.Begin(); itr != doc.End(); ++itr)
		{
			Polyhedron p;
			p.Deserialize(*itr);
			_polyhedrons.push_back(p);
		}

		return true;
	}

	bool Polyhedrons::Serialize(rapidjson::Writer<rapidjson::StringBuffer>* writer) const
	{
		writer->StartArray();

		for (std::vector<Polyhedron>::const_iterator it = _polyhedrons.begin(); it != _polyhedrons.end(); it++)
		{
			(*it).Serialize(writer);
		}
		writer->EndArray();

		return true;
	}
}