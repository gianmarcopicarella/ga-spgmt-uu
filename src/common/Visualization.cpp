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
		}

		void VisualizeLowerEnvelope(const LowerEnvelope3d& aLowerEnvelope, bool aShouldShowTrianglesFlag)
		{
			constexpr auto width = 1280.f;
			constexpr auto height = 720.f;
			constexpr auto vertexRadius = 3.f;

			std::vector<sf::Vector2f> positions;
			std::vector<sf::Vertex> edges;

			CGAL_precondition(!std::holds_alternative<std::monostate>(aLowerEnvelope));

			if (std::holds_alternative<std::vector<Edge<Point3>>>(aLowerEnvelope))
			{
				auto lowerEnvelope = std::get<std::vector<Edge<Point3>>>(aLowerEnvelope);

				if (lowerEnvelope.size() == 2 && lowerEnvelope[0].myType == EdgeType::LINE && lowerEnvelope[1].myType == EdgeType::LINE)
				{
					std::vector<Point3> points;
					points.emplace_back(lowerEnvelope[0].myStart);
					points.emplace_back(lowerEnvelope[0].myEnd);

					constexpr auto rangeOffset = 1.f;
					constexpr auto videoOffset = -50.f;
					const LinearMapToScreen mapping{ points, width, height, rangeOffset, videoOffset };
					Point2 cgalVertexPosition{ points.front().x(), points.front().y() };
					const auto& startVertexPosition = mapping.MapPointToScreen(cgalVertexPosition);
					cgalVertexPosition = Point2{ points.back().x(), points.back().y() };
					const auto& endVertexPosition = mapping.MapPointToScreen(cgalVertexPosition);
					edges.push_back(sf::Vertex{ startVertexPosition, sf::Color::White });
					edges.push_back(sf::Vertex{ endVertexPosition, sf::Color::White });
				}
				else
				{
					std::vector<Point3> uniquePoints;
					std::vector<Edge<Point3>> uniqueEdges;

					{
						const auto sortEdges = [](auto& aFirstEdge, auto& aSecondEdge)
						{
							return aFirstEdge.myStart < aSecondEdge.myStart ||
								(aFirstEdge.myStart == aSecondEdge.myStart && aFirstEdge.myEnd < aSecondEdge.myEnd);
						};

						const auto sameEdges = [](auto& aFirstEdge, auto& aSecondEdge)
						{
							return aFirstEdge.myStart == aSecondEdge.myStart && aFirstEdge.myEnd == aSecondEdge.myEnd;
						};

						for (const auto& edge : lowerEnvelope)
						{
							CGAL_precondition(edge.myType != EdgeType::LINE);

							if (edge.myType == EdgeType::SEGMENT ||
								edge.myType == EdgeType::HALF_EDGE_SF)
							{
								uniquePoints.emplace_back(edge.myStart);

								Edge<Point3> newEdge;
								newEdge.myType = edge.myType;
								newEdge.myLowestLeftPlane = edge.myLowestLeftPlane;
								newEdge.myStart = edge.myStart;
								newEdge.myEnd = edge.myEnd;

								if (newEdge.myType == EdgeType::SEGMENT)
								{
									if (newEdge.myStart > newEdge.myEnd)
									{
										std::swap(newEdge.myStart, newEdge.myEnd);
									}
									uniquePoints.emplace_back(edge.myEnd);
								}
								uniqueEdges.emplace_back(newEdge);
							}
						}

						if (aShouldShowTrianglesFlag)
						{
							for (const auto& edge : lowerEnvelope)
							{
								if (edge.myType == EdgeType::SEGMENT_TRIANGLE)
								{
									Edge<Point3> newEdge;
									newEdge.myType = edge.myType;
									newEdge.myLowestLeftPlane = edge.myLowestLeftPlane;
									newEdge.myStart = edge.myStart;
									newEdge.myEnd = edge.myEnd;
									if (newEdge.myStart > newEdge.myEnd)
									{
										std::swap(newEdge.myStart, newEdge.myEnd);
									}
									uniqueEdges.emplace_back(newEdge);
								}
							}
						}

						std::sort(uniquePoints.begin(), uniquePoints.end());
						std::sort(uniqueEdges.begin(), uniqueEdges.end(), sortEdges);
						uniquePoints.erase(std::unique(uniquePoints.begin(), uniquePoints.end()), uniquePoints.end());
						uniqueEdges.erase(std::unique(uniqueEdges.begin(), uniqueEdges.end(), sameEdges), uniqueEdges.end());
					}

					const LinearMapToScreen mapping{ uniquePoints, width, height };

					for (const auto& point : uniquePoints)
					{
						const auto& vertexPosition = mapping.MapPointToScreen(Point2{ point.x(), point.y() });
						positions.emplace_back(vertexPosition);
					}

					std::map<Point3, sf::Vector2f, CGAL::Less<Point3, Point3>> infinityMap;

					for (const auto& edge : uniqueEdges)
					{
						if (edge.myType < EdgeType::SEGMENT_TRIANGLE)
						{
							const auto& startPosition = mapping.MapPointToScreen(Point2{ edge.myStart.x(), edge.myStart.y() });
							edges.emplace_back(sf::Vertex{ startPosition });

							if (edge.myType == EdgeType::SEGMENT)
							{
								const auto& position = mapping.MapPointToScreen(Point2{ edge.myEnd.x(), edge.myEnd.y() });
								edges.emplace_back(sf::Vertex{ position });
							}
							else
							{
								const auto& infiniteVertexPosition =
									mapping.MapPointToScreen(Point2{ edge.myEnd.x(), edge.myEnd.y() });
								constexpr auto distance = 10000.f;
								const auto position = sf::Vector2f{ infiniteVertexPosition - startPosition }.normalized() * distance +
									startPosition;
								edges.emplace_back(sf::Vertex{ position });
							}

							if (aShouldShowTrianglesFlag)
							{
								infinityMap.insert(std::make_pair(edge.myEnd, edges[edges.size() - 1].position));
								infinityMap.insert(std::make_pair(edge.myStart, edges[edges.size() - 2].position));
							}
						}
					}

					if (aShouldShowTrianglesFlag)
					{
						for (const auto& edge : uniqueEdges)
						{
							if (edge.myType == EdgeType::SEGMENT_TRIANGLE)
							{
								edges.emplace_back(sf::Vertex{ infinityMap.at(edge.myStart), sf::Color::Red });
								edges.emplace_back(sf::Vertex{ infinityMap.at(edge.myEnd), sf::Color::Red });
							}
						}
					}
				}
			}
			else
			{
				CGAL_precondition(std::holds_alternative<size_t>(aLowerEnvelope));
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

				if (positions.size() > 0)
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