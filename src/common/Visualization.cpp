#include "Visualization.h"
#include "BruteForce.h"

#include <SFML/Graphics.hpp>

#include <CGAL/bounding_box.h>
#include <CGAL/assertions.h>

namespace SPGMT
{
	namespace Visualization
	{
		namespace
		{
			FT locLinearRangeMap(const FT& t, const FT& a, const FT& b, const FT& c, const FT& d)
			{
				return c + ((d - c) / (b - a)) * (t - a);
			}

			struct EdgePairHash
			{
				inline std::size_t operator()(const std::pair<int, int>& v) const
				{
					return v.first * 31 + v.second;
				}
			};
		}

		struct LinearMapToScreen
		{
			float myWidth, myHeight, myRangeOffset, myVideoOffset;
			Cube myBB;
			LinearMapToScreen(const std::vector<Point3>& somePoints, const float aWidth, const float anHeight, 
				const float aRangeOffset = 1.f, const float aVideoOffset = 50.f) :
				myWidth(aWidth), myHeight(anHeight), myRangeOffset(aRangeOffset), 
				myVideoOffset(aVideoOffset), myBB(CGAL::bounding_box(somePoints.begin(), somePoints.end()))
			{}
			sf::Vector2f MapPointToScreen(const Point2& aPosition) const
			{
				const auto x = locLinearRangeMap(aPosition.x(), myBB.xmin() - myRangeOffset, myBB.xmax() + myRangeOffset, 
					myVideoOffset, myWidth - myVideoOffset);
				const auto y = locLinearRangeMap(aPosition.y(), myBB.ymin() - myRangeOffset, myBB.ymax() + myRangeOffset, 
					myVideoOffset, myHeight - myVideoOffset);
				return sf::Vector2f{ (float)CGAL::to_double(x), myHeight - (float)CGAL::to_double(y) };
			}
		};

		void VisualizeLowerEnvelope(const std::vector<Vertex>& someVertices)
		{
			constexpr auto width = 1280.f;
			constexpr auto height = 720.f;
			constexpr auto vertexRadius = 3.f;

			std::vector<sf::Vector2f> positions;
			std::vector<sf::Vertex> edges;
			
			// If there are one or more edges with both ends at infinity
			// then handle this case properly
			constexpr auto singleEdgeFaceFind = [](const auto& aFace) { return aFace.myType == FaceType::UNBOUNDED_ONE_EDGE; };
			const auto& faces = ExtractLowerEnvelopeFaces(someVertices);
			const auto areThereSingleEdgeFaces = std::find_if(faces.begin(), faces.end(), singleEdgeFaceFind) != faces.end();
			if (areThereSingleEdgeFaces)
			{
				// Collect lines ends
				constexpr auto transform = [](const auto& aVertex) {return aVertex.myPoint; };
				std::vector<Point3> points;
				std::transform(someVertices.begin(), someVertices.end(), std::back_inserter(points), transform);
				// Build BB and map points to screen space
				constexpr auto rangeOffset = 1.f;
				constexpr auto videoOffset = -50.f;
				const LinearMapToScreen mapping{ points, width, height, rangeOffset, videoOffset };
				Point2 cgalVertexPosition{ someVertices.front().myPoint.x(), someVertices.front().myPoint.y()};
				const auto& startVertexPosition = mapping.MapPointToScreen(cgalVertexPosition);
				cgalVertexPosition = Point2{ someVertices.back().myPoint.x(), someVertices.back().myPoint.y()};
				const auto& endVertexPosition = mapping.MapPointToScreen(cgalVertexPosition);
				edges.push_back(sf::Vertex{ startVertexPosition, sf::Color::White });
				edges.push_back(sf::Vertex{ endVertexPosition, sf::Color::White });
			}
			else if(someVertices.size() > 1)
			{
				std::vector<Point3> points;
				{
					constexpr auto copy = [](const auto& aVertex) { return aVertex.myType == VertexType::FINITE; };
					constexpr auto transform = [](const auto& aVertex) {return aVertex.myPoint; };
					std::vector<Vertex> finiteVertices;
					std::copy_if(someVertices.begin(), someVertices.end(), std::back_inserter(finiteVertices), copy);
					std::transform(finiteVertices.begin(), finiteVertices.end(), std::back_inserter(points), transform);
				}
				const LinearMapToScreen mapping{ points, width, height};
				std::unordered_set<std::pair<int, int>, EdgePairHash> uniqueEdges;
				for (int i = 0; i < someVertices.size(); ++i)
				{
					if (someVertices[i].myType == VertexType::FINITE)
					{
						const auto& vertexPosition = 
							mapping.MapPointToScreen(Point2{ someVertices[i].myPoint.x(), someVertices[i].myPoint.y() });
						positions.push_back(vertexPosition);
						for (const auto k : someVertices[i].mySortedNeighboursIndices)
						{
							uniqueEdges.insert(std::make_pair(std::min(i, k), std::max(i, k)));
						}
					}
					else
					{
						CGAL_precondition(someVertices[i].mySortedNeighboursIndices.size() == 1);
						// Normalization is safe because the point at infinity is always far
						const auto neighbourIndex = someVertices[i].mySortedNeighboursIndices.front();
						const auto& finiteVertexPosition =
							mapping.MapPointToScreen(Point2{ someVertices[neighbourIndex].myPoint.x(), someVertices[neighbourIndex].myPoint.y() });
						const auto& infiniteVertexPosition =
							mapping.MapPointToScreen(Point2{ someVertices[i].myPoint.x(), someVertices[i].myPoint.y() });
						constexpr auto distance = 10000.f;
						positions.push_back(sf::Vector2f{ infiniteVertexPosition - finiteVertexPosition }.normalized() * distance + 
							finiteVertexPosition);
					}
				}
				for (const auto& edge : uniqueEdges)
				{
					edges.push_back(sf::Vertex{ positions[edge.first] });
					edges.push_back(sf::Vertex{ positions[edge.second] });
				}
			}
			else
			{
				CGAL_precondition(((someVertices.size() == 0) || (someVertices.size() == 1 && someVertices.front().myType == VertexType::INFINITE)));
			}

			sf::ContextSettings settings;
			settings.antialiasingLevel = 4;
			sf::View view(sf::Vector2f{ width, height } / 2.f, sf::Vector2f{ width, height });
			sf::RenderWindow window{ sf::VideoMode{ sf::Vector2u{(unsigned int)width, (unsigned int)height} }, "Lower Envelope 2D Visualization",
				sf::Style::Close | sf::Style::Titlebar, settings };
			window.setFramerateLimit(30);
			window.setView(view);

			const sf::Vector2f drawAdjust{ vertexRadius, vertexRadius };
			auto isLeftMousePressed{ false };
			sf::Vector2f mousePressOrigin;
			sf::Vector2i lastMousePos;
			sf::Vector2f lastViewCenter;
			sf::CircleShape circleShape(vertexRadius);
			circleShape.setFillColor(sf::Color::Magenta);

			while (window.isOpen())
			{
				for (auto event = sf::Event{}; window.pollEvent(event);)
				{
					switch (event.type)
					{
					case sf::Event::Closed:
						window.close();
						break;
					case sf::Event::MouseButtonPressed:
						if (event.mouseButton.button == sf::Mouse::Left)
						{
							isLeftMousePressed = true;
							mousePressOrigin.x = event.mouseButton.x;
							mousePressOrigin.y = event.mouseButton.y;
							lastViewCenter = view.getCenter();
						}
						break;
					case sf::Event::MouseButtonReleased:
						if (event.mouseButton.button == sf::Mouse::Left)
						{
							isLeftMousePressed = false;
						}
						break;
					case sf::Event::MouseWheelScrolled:
						if (event.mouseWheelScroll.wheel == sf::Mouse::VerticalWheel)
						{
							if (event.mouseWheelScroll.delta > 0)
							{
								view.zoom(0.97);
							}
							else
							{
								view.zoom(1.03);
							}
							window.setView(view);
						}
						break;
					default:
						break;
					}
				}
				window.clear();

				if (isLeftMousePressed)
				{
					sf::Vector2i mousePos = sf::Mouse::getPosition(window);
					if (mousePos != lastMousePos)
					{
						view.setCenter(lastViewCenter + mousePressOrigin - sf::Vector2f{ mousePos });
						window.setView(view);
						lastMousePos = mousePos;
					}
				}
				window.clear();

				if (edges.size() > 0)
				{
					window.draw(&edges[0], edges.size(), sf::PrimitiveType::Lines);
				}

				if (positions.size() > 2)
				{
					for (const auto& position : positions)
					{
						circleShape.setPosition(position - drawAdjust);
						window.draw(circleShape);
					}
				}

				window.display();
			}
		}
	}
}