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
		}

		struct LinearMapToScreen
		{
			FT myWidth, myHeight;
			FT myRangeOffset = 1.f;
			FT myVideoOffset = 50.f;
			Rec2 myBB;

			LinearMapToScreen(const std::vector<Vertex>& someVertices, const float aWidth, const float anHeight)
			{
				myWidth = aWidth;
				myHeight = anHeight;
				CGAL_precondition((myWidth != 0 && myHeight != 0 && myWidth > myVideoOffset && myHeight > myVideoOffset));
				CGAL_precondition(someVertices.size() > 0);
				std::vector<Point2> points;
				for (const auto& aVertex : someVertices)
				{
					if (aVertex.myType == VertexType::FINITE)
					{
						const Point2 point{ aVertex.myPoint.x(), aVertex.myPoint.y() };
						points.push_back(point);
					}
				}

				// Handling cases where there are less than 2 FINITE vertices in the Lower Envelope
				// Without these lines the bounding box function would trigger an assert
				if (points.size() == 1)
				{
					points.push_back(points.front());
				}

				if (points.size() > 0)
				{
					myBB = CGAL::bounding_box(points.begin(), points.end());
				}
			}

			sf::Vector2f MapPointToScreen(const Point2& aPosition) const
			{
				FT newX = locLinearRangeMap(aPosition.x(), myBB.xmin() - myRangeOffset, myBB.xmax() + myRangeOffset, myVideoOffset, myWidth - myVideoOffset);
				FT newY = locLinearRangeMap(aPosition.y(), myBB.ymin() - myRangeOffset, myBB.ymax() + myRangeOffset, myVideoOffset, myHeight - myVideoOffset);
				return sf::Vector2f{ (float)CGAL::to_double(newX), (float)CGAL::to_double(newY) };
			}

		};

		void VisualizeLowerEnvelope(const std::vector<Vertex>& someVertices)
		{
			std::vector<sf::Vertex> edges;
			std::vector<sf::Vector2f> verticesPositions;

			constexpr auto vertexRadius = 3.f;
			constexpr auto width = 1280.f;
			constexpr auto height = 720.f;
			
			if (someVertices.size() > 1)
			{
				struct pairHash
				{
					inline std::size_t operator()(const std::pair<int, int>& v) const
					{
						return v.first * 31 + v.second;
					}
				};

				std::unordered_set<std::pair<int, int>, pairHash> uniqueEdges;
				std::vector<sf::Vector2f> sfmlData;
				const LinearMapToScreen mapping{ someVertices, width, height };
				const sf::Vector2f drawAdjust{ vertexRadius, vertexRadius };

				for (int i = 0; i < someVertices.size(); ++i)
				{
					if (someVertices[i].myType == VertexType::FINITE)
					{
						sfmlData.push_back(
							mapping.MapPointToScreen(Point2{ someVertices[i].myPoint.x(), someVertices[i].myPoint.y() }));
						for (const auto k : someVertices[i].mySortedNeighboursIndices)
						{
							uniqueEdges.insert(std::make_pair(std::min(i, k), std::max(i, k)));
						}
						verticesPositions.push_back(sfmlData.back() - drawAdjust);
					}
					else
					{
						CGAL_precondition(someVertices[i].mySortedNeighboursIndices.size() == 1);
						// Normalization is safe because the point at infinity is always far
						const auto neighbourIndex = someVertices[i].mySortedNeighboursIndices.front();
						const Vec3 halfEdgeDir{ someVertices[neighbourIndex].myPoint, someVertices[i].myPoint };
						sfmlData.push_back(sf::Vector2f{
							(float)CGAL::to_double(halfEdgeDir.x()),
							(float)CGAL::to_double(halfEdgeDir.y())
							}.normalized());
						uniqueEdges.insert(std::make_pair(std::min(i, neighbourIndex), std::max(i, neighbourIndex)));
					}
				}

				for (const auto& edge : uniqueEdges)
				{
					if (someVertices[edge.first].myType == VertexType::FINITE &&
						someVertices[edge.second].myType == VertexType::FINITE)
					{
						edges.push_back(sf::Vertex{ sfmlData[edge.first] });
						edges.push_back(sf::Vertex{ sfmlData[edge.second] });
					}
					else if (someVertices[edge.first].myType == VertexType::FINITE)
					{
						edges.push_back(sf::Vertex{ sfmlData[edge.first] });
						edges.push_back(sf::Vertex{ sfmlData[edge.first] + sfmlData[edge.second] * 10000.f });
					}
					else if (someVertices[edge.second].myType == VertexType::FINITE)
					{
						edges.push_back(sf::Vertex{ sfmlData[edge.second] });
						edges.push_back(sf::Vertex{ sfmlData[edge.second] + sfmlData[edge.first] * 10000.f });
					}
					else
					{
						const sf::Vector2f origin{ 640,360 };
						edges.push_back(sf::Vertex{ origin + sfmlData[edge.second] * 10000.f });
						edges.push_back(sf::Vertex{ origin + sfmlData[edge.first] * 10000.f });
					}
				}
			}
			else
			{
				CGAL_precondition((someVertices.size() == 1 && someVertices.front().myType == VertexType::INFINITE));
			}

			sf::ContextSettings settings;
			settings.antialiasingLevel = 4;

			sf::View view(sf::Vector2f{ width, height } / 2.f, sf::Vector2f{ width, height });
			sf::RenderWindow window{ sf::VideoMode{ sf::Vector2u{(unsigned int)width, (unsigned int)height} }, "Lower Envelope 2D Visualization",
				sf::Style::Close | sf::Style::Titlebar, settings };
			window.setFramerateLimit(30);
			window.setView(view);

			auto isLeftMousePressed{ false };
			auto isMouseWheelMoving{ false };
			sf::Vector2f mousePressOrigin;
			sf::Vector2i lastMousePos;
			sf::Vector2f lastViewCenter;

			sf::CircleShape circleShape(vertexRadius);
			circleShape.setFillColor(sf::Color::Magenta);
			//circleShape.setPointCount(100);

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
								view.zoom(0.98);
							}
							else
							{
								view.zoom(1.02);
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

				if (edges.size() > 0)
				{
					window.draw(&edges[0], edges.size(), sf::PrimitiveType::Lines);
				}

				if (verticesPositions.size() > 0)
				{
					for (const auto& position : verticesPositions) 
					{
						circleShape.setPosition(position);
						window.draw(circleShape);
					}
				}

				window.display();
			}
		}
	}
}