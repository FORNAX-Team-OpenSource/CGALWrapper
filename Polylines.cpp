#include <string>
#include "Polylines.h"

namespace JSONModels
{
	Polylines::Polylines()
	{}

	Polylines::~Polylines()
	{}

	bool Polylines::Deserialize(const rapidjson::Value& obj)
	{
		return true;
	}

	bool Polylines::Serialize(rapidjson::Writer<rapidjson::StringBuffer>* writer) const
	{
		writer->StartArray();
		for (std::vector<Polyline>::const_iterator it = _polylines.begin(); it != _polylines.end(); it++)
		{
			(*it).Serialize(writer);
		}
		writer->EndArray();
		return true;
	}
}
