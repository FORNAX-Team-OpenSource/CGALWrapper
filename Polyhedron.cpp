#include <string>
#include "Polyhedron.h"

namespace JSONModels
{
	Polyhedron::Polyhedron()
	{}

	Polyhedron::~Polyhedron()
	{}



	bool Polyhedron::Deserialize(const rapidjson::Value& obj)
	{
		//Area(obj["Area"].GetDouble());

		for (rapidjson::Value::ConstValueIterator itr = obj["Faces"].Begin(); itr != obj["Faces"].End(); ++itr)
		{
			Face p;
			p.Deserialize(*itr);
			_faces.push_back(p);
		}

		return true;
	}

	bool Polyhedron::Serialize(rapidjson::Writer<rapidjson::StringBuffer>* writer) const
	{
		writer->StartObject();

		writer->String("Faces");
		writer->StartArray();

		for (std::vector<Face>::const_iterator it = _faces.begin(); it != _faces.end(); it++)
		{
			(*it).Serialize(writer);
		}
		writer->EndArray();

		//writer->String("Area");
		//writer->Double(_area);

		writer->EndObject();

		return true;
	}
}
