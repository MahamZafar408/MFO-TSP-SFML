#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <limits>
#include <random>
#include <fstream>

using namespace std;


class FitnessFunction {
public:
    static constexpr double R = 6371.0; 

    static double haversine(double lat1, double lon1, double lat2, double lon2) {
        const double R = 6371.0; 
        double dLat = (lat2 - lat1) * 3.14159265358979323846 / 180.0;
        double dLon = (lon2 - lon1) * 3.14159265358979323846 / 180.0;
        lat1 = lat1 * 3.14159265358979323846 / 180.0;
        lat2 = lat2 * 3.14159265358979323846 / 180.0;

        double a = sin(dLat / 2) * sin(dLat / 2) +
            sin(dLon / 2) * sin(dLon / 2) * cos(lat1) * cos(lat2);
        double c = 2 * atan2(sqrt(a), sqrt(1 - a));
        return R * c;
    }

    static double computePathDistance(const vector<int>& path, const vector<vector<double>>& distances) {
        double totalDistance = 0;
        for (size_t i = 1; i < path.size(); ++i) {
            totalDistance += distances[path[i - 1]][path[i]];
        }
        totalDistance += distances[path.back()][path.front()];
        return totalDistance;
    }
};



class MothFlameOptimizer {
private:
    const vector<vector<double>>& distances;
    int populationSize;
    int maxIterations;
    vector<vector<double>>& allIterationDistances;
    vector<int> bestPath;
    double bestDistance;

public:
    MothFlameOptimizer(const vector<vector<double>>& distances, int populationSize, int maxIterations,
        vector<vector<double>>& allIterationDistances)
        : distances(distances), populationSize(populationSize), maxIterations(maxIterations),
        allIterationDistances(allIterationDistances), bestDistance(numeric_limits<double>::infinity()) {}

    vector<int> optimize(sf::RenderWindow& window, const vector<sf::Vector2f>& positions, sf::Font& font) {
        int n = distances.size();
        vector<vector<int>> moths(populationSize, vector<int>(n));
        random_device rd;
        mt19937 gen(42); 
        uniform_int_distribution<> dis(0, n - 1);

        
        for (auto& moth : moths) {
            iota(moth.begin(), moth.end(), 0);
            shuffle(moth.begin(), moth.end(), gen);
        }


        bestPath = moths[0];
        bestDistance = FitnessFunction::computePathDistance(bestPath, distances);

        ofstream csvFile("iteration_distances.csv");
        csvFile << "Iteration";
        for (int i = 0; i < populationSize; ++i) {
            csvFile << ",Moth " << (i + 1) << " Distance";
        }
        csvFile << "\n";

        for (int iteration = 0; iteration < maxIterations; ++iteration) {
            vector<double> distancesThisIteration;

            
            for (auto& moth : moths) {
                double distance = FitnessFunction::computePathDistance(moth, distances);
                distancesThisIteration.push_back(distance);
                if (distance < bestDistance) {
                    bestDistance = distance;
                    bestPath = moth;
                }
            }

            allIterationDistances.push_back(distancesThisIteration);

            csvFile << iteration + 1;
            for (double d : distancesThisIteration) {
                csvFile << "," << d;
            }
            csvFile << "\n";

            
               double t = (double)iteration / maxIterations;
               double b = 1.5;
            for (size_t i = 0; i < moths.size(); ++i) {
                if (moths[i] != bestPath) {
                    
                    vector<int> flame = bestPath;

                    double d = FitnessFunction::computePathDistance(moths[i], distances) -
                        FitnessFunction::computePathDistance(flame, distances);

                    for (int j = 0; j < n; ++j) {
                       
                        double move = b * d * cos(2 * 3.14159265358979323846 * t); 
                        int pos = (j + (int)(move * n)) % n; 

                        if (pos < 0) pos += n; 
                        swap(moths[i][j], moths[i][pos]);

                    }
                    

                   
                }
            }
            
        }

        csvFile.close();
        return bestPath;
    }

    double getBestDistance() const { return bestDistance; }
    const vector<int>& getBestPath() const { return bestPath; }
};

enum DisplayState {
    ADJACENCY_MATRIX,
    CITY_GRAPH,
    ITERATION_GRAPH,
    BEST_PATH,
    END
};
const vector<pair<string, pair<double, double>>> locations = {
    {"Islamabad", {33.6844, 73.0479}},
    {"Karachi", {24.8607, 67.0011}},
    {"Lahore", {31.5497, 74.3436}},
    {"Peshawar", {34.0151, 71.5249}},
    {"Quetta", {30.1798, 66.9750}},
    {"Multan", {30.1575, 71.5249}}
};
void drawAdjacencyMatrix(sf::RenderWindow& window, const vector<vector<double>>& distances, sf::Font& font) {
    const int margin = 20;
    const int cellSize = 100; 
    const int offsetX = 500; 
    const int offsetY = 250; 

    for (size_t i = 0; i <= distances.size(); ++i) {
        for (size_t j = 0; j <= distances.size(); ++j) {
            if (i == 0 && j == 0) continue; 

            sf::RectangleShape cell(sf::Vector2f(cellSize, cellSize));
            cell.setFillColor(sf::Color::White);
            cell.setOutlineThickness(1);
            cell.setOutlineColor(sf::Color::Black);
            cell.setPosition(offsetX + j * cellSize, offsetY + i * cellSize);
            window.draw(cell);

            sf::Text label;
            if (i == 0) {
                label.setString(locations[j - 1].first);
            }
            else if (j == 0) {
                label.setString(locations[i - 1].first);
            }
            else {
                stringstream ss;
                ss << fixed << setprecision(1) << distances[i - 1][j - 1];
                label.setString(ss.str());
            }

            label.setFont(font);
            label.setCharacterSize(18);
            label.setFillColor(sf::Color::Black);
            label.setPosition(offsetX + j * cellSize + 10, offsetY + i * cellSize + 10);
            window.draw(label);
        }
    }
}

void drawCityGraph(sf::RenderWindow& window, const vector<sf::Vector2f>& positions, const vector<vector<double>>& distances, sf::Font& font) {
    int n = positions.size();

    sf::Vector2u windowSize = window.getSize();
    float centerX = windowSize.x / 2.0f;
    float centerY = windowSize.y / 2.0f;

    float layoutRadius = min(windowSize.x, windowSize.y) / 3.0f;

    vector<sf::Vector2f> scaledPositions(n);
    for (int i = 0; i < n; ++i) {
        float angle = i * (2 * 3.14159265358979323846 / n);
        scaledPositions[i].x = centerX + cos(angle) * layoutRadius;
        scaledPositions[i].y = centerY + sin(angle) * layoutRadius;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            sf::Vertex line[] = {
                sf::Vertex(scaledPositions[i], sf::Color::Black),
                sf::Vertex(scaledPositions[j], sf::Color::Black)
            };
            window.draw(line, 2, sf::Lines);

            sf::Text distanceText;
            stringstream ss;
            ss << fixed << setprecision(1) << distances[i][j];
            distanceText.setString(ss.str() + " km");
            distanceText.setFont(font);
            distanceText.setCharacterSize(16);
            distanceText.setFillColor(sf::Color::Black);
            distanceText.setPosition(
                (scaledPositions[i].x + scaledPositions[j].x) / 2,
                (scaledPositions[i].y + scaledPositions[j].y) / 2
            );
            window.draw(distanceText);
        }
    }

    const int radius = 25;
    for (size_t i = 0; i < n; ++i) {
        sf::CircleShape node(radius);
        node.setFillColor(sf::Color::Blue);
        node.setPosition(scaledPositions[i].x - radius, scaledPositions[i].y - radius);
        window.draw(node);

        sf::Text cityName;
        cityName.setString(locations[i].first);
        cityName.setFont(font);
        cityName.setCharacterSize(18);
        cityName.setFillColor(sf::Color::Black);
        cityName.setPosition(scaledPositions[i].x - radius, scaledPositions[i].y - radius - 30);
        window.draw(cityName);
    }
}

void drawAllIterationGraph(sf::RenderWindow& window, const vector<vector<double>>& allIterationDistances, sf::Font& font) {
    
    const float windowWidth = window.getSize().x;
    const float windowHeight = window.getSize().y;

    const float graphWidth = 900.0f;  
    const float graphHeight = 500.0f; 
    const float marginLeft = (windowWidth - graphWidth) / 2;
    const float marginTop = (windowHeight - graphHeight) / 2; 

    double minDistance = numeric_limits<double>::max();
    double maxDistance = 0;

    for (const auto& iteration : allIterationDistances) {
        for (double d : iteration) {
            if (d < minDistance) minDistance = d;
            if (d > maxDistance) maxDistance = d;
        }
    }

    sf::Vertex xAxis[] = {
        sf::Vertex(sf::Vector2f(marginLeft, marginTop + graphHeight), sf::Color::Black),
        sf::Vertex(sf::Vector2f(marginLeft + graphWidth, marginTop + graphHeight), sf::Color::Black)
    };

    sf::Vertex yAxis[] = {
        sf::Vertex(sf::Vector2f(marginLeft, marginTop + graphHeight), sf::Color::Black),
        sf::Vertex(sf::Vector2f(marginLeft, marginTop), sf::Color::Black)
    };

    window.draw(xAxis, 2, sf::Lines);
    window.draw(yAxis, 2, sf::Lines);

    int xTicks = 10;
    int yTicks = 10; 

 
    for (int i = 0; i <= xTicks; ++i) {
        float xPos = marginLeft + (i * graphWidth / xTicks);
        sf::Vertex tick[] = {
            sf::Vertex(sf::Vector2f(xPos, marginTop + graphHeight), sf::Color::Black),
            sf::Vertex(sf::Vector2f(xPos, marginTop + graphHeight + 5), sf::Color::Black)
        };
        window.draw(tick, 2, sf::Lines);

        sf::Text label;
        label.setFont(font);
        label.setCharacterSize(14);
        label.setFillColor(sf::Color::Black);
        label.setString(to_string(i * (allIterationDistances.size() / xTicks)));
        label.setPosition(xPos - 10, marginTop + graphHeight + 10);
        window.draw(label);
    }

    for (int i = 0; i <= yTicks; ++i) {
        float yPos = marginTop + graphHeight - (i * graphHeight / yTicks);
        sf::Vertex tick[] = {
            sf::Vertex(sf::Vector2f(marginLeft - 5, yPos), sf::Color::Black),
            sf::Vertex(sf::Vector2f(marginLeft, yPos), sf::Color::Black)
        };
        window.draw(tick, 2, sf::Lines);

        double distanceValue = minDistance + i * (maxDistance - minDistance) / yTicks;
        sf::Text label;
        label.setFont(font);
        label.setCharacterSize(14);
        label.setFillColor(sf::Color::Black);
        stringstream ss;
        ss << fixed << setprecision(2) << distanceValue;
        label.setString(ss.str());
        label.setPosition(marginLeft - 50, yPos - 10);
        window.draw(label);
    }

    sf::VertexArray lineGraph(sf::LinesStrip);

    for (size_t i = 0; i < allIterationDistances.size(); ++i) {
        double avgDistance = accumulate(allIterationDistances[i].begin(), allIterationDistances[i].end(), 0.0) / allIterationDistances[i].size();

        float xPos = marginLeft + (i * graphWidth / allIterationDistances.size());
        float yPos = marginTop + graphHeight - ((avgDistance - minDistance) / (maxDistance - minDistance) * graphHeight);

        lineGraph.append(sf::Vertex(sf::Vector2f(xPos, yPos), sf::Color::Red));
    }

    window.draw(lineGraph);
}


void waitForEnter(sf::RenderWindow& window) {
    while (true) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Enter) {
                return;
            }
        }
    }
}

int main() {

    const int windowWidth = 800;
    const int windowHeight = 600;

    int n = locations.size();
    vector<vector<double>> distances(n, vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            distances[i][j] = FitnessFunction::haversine(
                locations[i].second.first, locations[i].second.second,
                locations[j].second.first, locations[j].second.second);
        }
    }

    
    vector<sf::Vector2f> positions;
    for (int i = 0; i < n; ++i) {
        float angle = i * (2 * 3.14159265358979323846 / n);
        float x = windowWidth / 2 + cos(angle) * (windowWidth / 3);
        float y = windowHeight / 2 + sin(angle) * (windowHeight / 3);
        positions.emplace_back(x, y);
    }

    sf::RenderWindow window(sf::VideoMode::getDesktopMode(), "City Adjacency Graph", sf::Style::Default);
    sf::Font font;
    if (!font.loadFromFile("arial.ttf")) {
        cerr << "Failed to load font.\n";
        return -1;
    }
    vector<vector<double>> allIterationDistances;
    MothFlameOptimizer mfo(distances, 18, 2000, allIterationDistances);
    vector<int> bestPath = mfo.optimize(window, positions, font);
    double bestDistance = mfo.getBestDistance();

    DisplayState currentState = ADJACENCY_MATRIX;

    while (window.isOpen()) {
        window.clear(sf::Color::White);

        switch (currentState) {
        case ADJACENCY_MATRIX:
            drawAdjacencyMatrix(window, distances, font);
            break;
        case CITY_GRAPH:
            drawCityGraph(window, positions, distances, font);
            break;
        case ITERATION_GRAPH:
            drawAllIterationGraph(window, allIterationDistances, font);
            break;
        case BEST_PATH: {
            const int radius = 20;

            for (size_t i = 0; i < bestPath.size(); ++i) {
                int from = bestPath[i];
                int to = bestPath[(i + 1) % bestPath.size()]; // Loop back to the first city

                sf::Vertex line[] = {
                    sf::Vertex(positions[from], sf::Color::Black),
                    sf::Vertex(positions[to], sf::Color::Black)
                };
                window.draw(line, 2, sf::Lines);

                sf::Text distanceText;
                stringstream ss;
                ss << fixed << setprecision(1) << distances[from][to] << " km";
                distanceText.setString(ss.str());
                distanceText.setFont(font);
                distanceText.setCharacterSize(14);
                distanceText.setFillColor(sf::Color::Black);
                distanceText.setPosition(
                    (positions[from].x + positions[to].x) / 2,
                    (positions[from].y + positions[to].y) / 2
                );
                window.draw(distanceText);
            }
            for (size_t i = 0; i < bestPath.size(); ++i) {
                int cityIndex = bestPath[i];
                sf::CircleShape node(radius);
                node.setFillColor(sf::Color::Blue);
                node.setPosition(positions[cityIndex].x - radius, positions[cityIndex].y - radius);
                window.draw(node);

                sf::Text cityName;
                cityName.setString(locations[cityIndex].first);
                cityName.setFont(font);
                cityName.setCharacterSize(14);
                cityName.setFillColor(sf::Color::Black);
                cityName.setPosition(positions[cityIndex].x - radius, positions[cityIndex].y - radius - 20);
                window.draw(cityName);
            }

            sf::Text distanceText, pathText;
            distanceText.setString("Best Distance: " + to_string(bestDistance) + " km");
            distanceText.setFont(font);
            distanceText.setCharacterSize(16);
            distanceText.setFillColor(sf::Color::Black);
            distanceText.setPosition(50, 500);
            window.draw(distanceText);

            stringstream pathStream;
            for (size_t i = 0; i < bestPath.size(); ++i) {
                pathStream << locations[bestPath[i]].first;
                if (i != bestPath.size() - 1) {
                    pathStream << " -> ";
                }
            }
            pathStream << " -> " << locations[bestPath[0]].first; // Close the loop

            pathText.setString("Best Path: " + pathStream.str());
            pathText.setFont(font);
            pathText.setCharacterSize(16);
            pathText.setFillColor(sf::Color::Black);
            pathText.setPosition(50, 530);
            window.draw(pathText);
            break;
        }
        case END:
            window.close();
            break;
        }

        window.display();
        waitForEnter(window);

        if (currentState == END) {
            window.close();
        }
        else {
            currentState = static_cast<DisplayState>(static_cast<int>(currentState) + 1);
        }
    }

    return 0;
}

