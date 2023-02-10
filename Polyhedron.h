#pragma once

#include "JSONBase.h"
#include <list>
#include <vector>
#include <string>
#include "Face.h"

namespace JSONModels
{
	class Polyhedron : public JSONBase
	{
	public:
		Polyhedron();
		virtual ~Polyhedron();

		virtual bool Deserialize(const rapidjson::Value& obj);
		virtual bool Serialize(rapidjson::Writer<rapidjson::StringBuffer>* writer) const;

		// Getters/Setters.

		//double Area() const { return _area; }
		//void Area(double area) { _area = area; }

		std::vector<Face>& FaceList() { return _faces; }
		void FaceList(std::vector<Face> faces) { _faces = faces; }
		/*VerticesWithHoles verticesWithHoles() const { return _verticesWithHoles; }
		void verticesWithHoles(VerticesWithHoles verticesWithHoles) { _verticesWithHoles = verticesWithHoles; }*/

		//

	private:
		//double _area;
		std::vector<Face> _faces;
		//VerticesWithHoles _verticesWithHoles;

	};
}

