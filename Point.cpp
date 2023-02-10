#include <string>
#include "Point.h"

namespace JSONModels
{
	Point::Point()
	{}

	Point::~Point()
	{}



	bool Point::Deserialize(const rapidjson::Value& obj)
	{
		X(obj["X"].GetDouble());
		Y(obj["Y"].GetDouble());
		Z(obj["Z"].GetDouble());

		return true;
	}

	bool Point::Serialize(rapidjson::Writer<rapidjson::StringBuffer>* writer) const
	{
		writer->StartObject();

		writer->String("X");
		writer->Double(_x);

		writer->String("Y");
		writer->Double(_y);

		writer->String("Z");
		writer->Double(_z);

		writer->EndObject();

		return true;
	}
}
