#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <CanvasPoint.h>
#include <Colour.h>
#include <CanvasTriangle.h>
#include <TextureMap.h>
#include <TexturePoint.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <unordered_map>
#define WIDTH 500
#define HEIGHT 400
glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 4.0);
glm::mat3 cameraOrientation = glm::mat3(1, 0 , 0, 0, 1, 0, 0 , 0, 1);
glm::vec3 lightSource = glm::vec3(0, 0.35, 0);
bool orbit = false;
int renderType = 1;

 
void drawLine(DrawingWindow &window, CanvasPoint to, CanvasPoint from, uint32_t colour) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float steps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff/steps;
	float yStepSize = yDiff/steps;
	for (float i=0.0; i<steps; i++) {
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);
		if(x >= 0 && round(x) <WIDTH){
			if(y >= 0 && round(y) < HEIGHT){
				window.setPixelColour(round(x), round(y), colour);
			}
		}
	}
}

void drawObjectLine(DrawingWindow &window, CanvasPoint to, CanvasPoint from, uint32_t colour, float depthArray[WIDTH][HEIGHT]) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float depthDiff = to.depth - from.depth;
	float steps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff/steps;
	float yStepSize = yDiff/steps;
	float depthStepSize = depthDiff/steps;
	for (float i=0.0; i<steps; i++) {
		if(i == 0) {
			float xf = floor(from.x + (xStepSize*i));
			float yf= floor(from.y + (yStepSize*i));
			int x = static_cast<int>(xf);
			int y = static_cast<int>(yf);
			if(x >= 0 && x <WIDTH){
				if(y >= 0 && y < HEIGHT){
					float depth  = -1/(from.depth + (depthStepSize*i));
					if (depth > depthArray[x][y]){
						window.setPixelColour(x, y, colour);
						depthArray[x][y] = depth;
					}
				}
			}
		}else if(i == (steps - 1)) {
			float xf = ceil(from.x + (xStepSize*i));
			float yf= ceil(from.y + (yStepSize*i));
			int x = static_cast<int>(xf);
			int y = static_cast<int>(yf);
			if(x >= 0 && x <WIDTH){
				if(y >= 0 && y < HEIGHT){
					float depth  = -1/(from.depth + (depthStepSize*i));
					if (depth > depthArray[x][y]) {
						window.setPixelColour(x, y, colour);
						depthArray[x][y] = depth;
					}
				}
			}
		}else {
			float xf = round(from.x + (xStepSize*i));
			float yf= round(from.y + (yStepSize*i));
			int x = static_cast<int>(xf);
			int y = static_cast<int>(yf);
			if(x >= 0 && x <WIDTH){
				if(y >= 0 && y < HEIGHT){
					float depth  = -1/(from.depth + (depthStepSize*i));
					if (depth > depthArray[x][y]) {
						window.setPixelColour(x, y, colour);
						depthArray[x][y] = depth;
					}
				}
			}
		}
		
	}
}

std::vector<CanvasPoint> interpolateLine(CanvasPoint to, CanvasPoint from) {
	std::vector<CanvasPoint> line;
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float depthDiff = to.depth - from.depth;
	float steps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff/steps;
	float yStepSize = yDiff/steps;
	float depthStepSize = depthDiff/steps;
	for (float i=0.0; i<steps; i++) {
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);
		float depth = from.depth + (depthStepSize*i);
		CanvasPoint point = CanvasPoint(round(x), round(y));
		point.depth = depth;
		line.push_back(point);
	}
	return line;
}

std::vector<CanvasPoint> interpolateLineTex(TexturePoint to, TexturePoint from) {
	std::vector<CanvasPoint> line;
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float steps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff/steps;
	float yStepSize = yDiff/steps;
	for (float i=0.0; i<steps; i++) {
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);
		line.push_back(CanvasPoint(round(x),round(y)));
	}
	return line;
}

std::vector<CanvasPoint> interpolateTextureLine(CanvasPoint to, CanvasPoint from, float steps) {
	std::vector<CanvasPoint> line;
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float xStepSize = xDiff/steps;
	float yStepSize = yDiff/steps;
	for (float i=0.0; i<steps; i++) {
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);
		line.push_back(CanvasPoint(round(x),round(y)));
	}
	return line;
}

void drawTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colourGiven){
	uint32_t colour = (255 << 24) + (colourGiven.red << 16) + (colourGiven.green << 8) + colourGiven.blue;
	drawLine(window, triangle.v0(), triangle.v1(), colour);
	drawLine(window, triangle.v0(), triangle.v2(), colour);
	drawLine(window, triangle.v2(), triangle.v1(), colour);
}

std::vector<CanvasPoint> sortVertices(CanvasTriangle triangle){
	CanvasPoint v0 = triangle.v0();
	CanvasPoint v1 = triangle.v1();
	CanvasPoint v2 = triangle.v2();
	std::vector<CanvasPoint> vertices;
	if (v1.y < v0.y){ 
		std::swap(v0, v1);
		if (v2.y < v1.y) {
			std::swap(v2, v1);
			if(v1.y < v0.y) std::swap(v1, v0);
		}
	}
	else if (v2.y < v1.y) {
		std::swap(v2, v1);
		if(v1.y < v0.y) std::swap(v1, v0);
	}
	vertices.push_back(v0);
	vertices.push_back(v1);
	vertices.push_back(v2);
	return vertices;
}

CanvasPoint calculateMidpoint(std::vector<CanvasPoint> vertices) {
	float midX;
	float midY;
	midY = vertices[1].y;
	float ratio = (vertices[1].y - vertices[0].y)/(vertices[2].y - vertices[0].y);
	midX = vertices[0].x + (ratio * (vertices[2].x - vertices[0].x));
	float midDepth = vertices[0].depth + (ratio * (vertices[2].depth - vertices[0].depth));                                                                                                                                                
	CanvasPoint point = CanvasPoint(round(midX), midY);
	point.depth = midDepth;
	return point;
}

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colourGiven) {
	uint32_t colour = (255 << 24) + (colourGiven.red << 16) + (colourGiven.green << 8) + colourGiven.blue;
	std::vector<CanvasPoint> sortedVertices = sortVertices(triangle);
	CanvasPoint midPoint = calculateMidpoint(sortedVertices);  
	std::vector<CanvasPoint> topLine1 = interpolateLine(sortedVertices[0], midPoint);
	std::vector<CanvasPoint> topLine2 = interpolateLine(sortedVertices[0], sortedVertices[1]);
	for(uint32_t i = 0; i < topLine1.size(); i++) {
		for(uint32_t j = 0; j < topLine2.size(); j++) {
			if (topLine2[j].y == topLine1[i].y){
				drawLine(window, topLine1[i], topLine2[j], colour);
			}
		}
	}
	std::vector<CanvasPoint> bottomLine1 = interpolateLine(sortedVertices[2], midPoint);
	std::vector<CanvasPoint> bottomLine2 = interpolateLine(sortedVertices[2], sortedVertices[1]);
	for(uint32_t i = 0; i < bottomLine1.size(); i++) {
		for(uint32_t j = 0; j < bottomLine2.size(); j++) {
			if (bottomLine2[j].y == bottomLine1[i].y){
				drawLine(window, bottomLine1[i], bottomLine2[j], colour);
			}
		}
	}
}

CanvasPoint calculatetexturedMidpoint(std::vector<CanvasPoint> vertices) {
	float midX;
	float midY;
	float texMidX;
	float texMidY;
	midY = vertices[1].y;
	float ratio = (vertices[1].y - vertices[0].y)/(vertices[2].y - vertices[0].y);
	midX = vertices[0].x + (ratio * (vertices[2].x - vertices[0].x));
	texMidX = vertices[0].texturePoint.x + (ratio * (vertices[2].texturePoint.x - vertices[0].texturePoint.x));
	texMidY = vertices[0].texturePoint.y + (ratio * (vertices[2].texturePoint.y - vertices[0].texturePoint.y));
	CanvasPoint midPoint = CanvasPoint(round(midX), midY);
	midPoint.texturePoint = TexturePoint(round(texMidX), round(texMidY)); 
	return midPoint; 
}

void drawTexturedLine(DrawingWindow &window, CanvasPoint to, CanvasPoint from, CanvasPoint texTo, CanvasPoint texFrom, TextureMap texture) {
	float xDiff = to.x - from.x;
	float steps = xDiff;
	float xStepSize = 1;
	std::vector<CanvasPoint> textureLine = interpolateTextureLine(texTo, texFrom, steps);
	for (float i=0.0; i<steps; i++) {
		float x = from.x + (xStepSize*i);
		float y = from.y;
		float pixel = (textureLine[i].y * texture.width) + textureLine[i].x;
		uint32_t colour = texture.pixels[pixel];
		window.setPixelColour(round(x), round(y), colour);
	}
}

void drawTexturedTriangle(DrawingWindow &window, std::vector<CanvasPoint> triangle, TextureMap texture, std::vector<TexturePoint> texturePoints) {
	triangle[0].texturePoint = texturePoints[0];
	triangle[1].texturePoint = texturePoints[1];
	triangle[2].texturePoint = texturePoints[2];
	CanvasTriangle triangle1 = CanvasTriangle(triangle[0], triangle[1], triangle[2]);
	std::vector<CanvasPoint> sortedVertices = sortVertices(triangle1);
	CanvasPoint midPoint = calculatetexturedMidpoint(sortedVertices);
	std::vector<CanvasPoint> topLine1 = interpolateLine(sortedVertices[0], midPoint);
	std::vector<CanvasPoint> topLine2 = interpolateLine(sortedVertices[0], sortedVertices[1]);
	std::vector<CanvasPoint> texTopLine1 = interpolateLineTex(sortedVertices[0].texturePoint, midPoint.texturePoint);
	std::vector<CanvasPoint> texTopLine2 = interpolateLineTex(sortedVertices[0].texturePoint, sortedVertices[1].texturePoint);
	for(float i = 0; i < topLine1.size(); i++) {
		for(float j = 0; j < topLine2.size(); j++) {
			if (topLine2[j].y == topLine1[i].y){
				float ratio1 = i/topLine1.size();
				float ratio2 = j/topLine2.size();
				float k = round(texTopLine1.size()*ratio1);
				float l = round(texTopLine2.size()*ratio2);
				drawTexturedLine(window, topLine1[i], topLine2[j], texTopLine1[k], texTopLine2[l], texture);
			}
		}
	}
	std::vector<CanvasPoint> bottomLine1 = interpolateLine(midPoint, sortedVertices[2]);
	std::vector<CanvasPoint> bottomLine2 = interpolateLine(sortedVertices[1], sortedVertices[2]);
	std::vector<CanvasPoint> texBottomLine1 = interpolateLineTex(midPoint.texturePoint, sortedVertices[2].texturePoint);
	std::vector<CanvasPoint> texBottomLine2 = interpolateLineTex(sortedVertices[1].texturePoint, sortedVertices[2].texturePoint);

	for(float i = 0; i < bottomLine1.size(); i++) {
		for(float j = 0; j < bottomLine2.size(); j++) {
			if (bottomLine2[j].y == bottomLine1[i].y){
				float ratio1 = i/bottomLine1.size();
				float ratio2 = j/bottomLine2.size();
				float k = round(texBottomLine1.size()*ratio1);
				float l = round(texBottomLine2.size()*ratio2);
				drawTexturedLine(window, bottomLine1[i], bottomLine2[j], texBottomLine1[k], texBottomLine2[l], texture);
			}
		}
	}
	drawTriangle(window, triangle1, Colour(255,255,255));
}

CanvasTriangle randTriangle() {
	CanvasPoint v1 = CanvasPoint(rand()%WIDTH, rand()%HEIGHT);
	CanvasPoint v2 = CanvasPoint(rand()%WIDTH, rand()%HEIGHT);
	CanvasPoint v3 = CanvasPoint(rand()%WIDTH, rand()%HEIGHT);
	CanvasTriangle triangle = CanvasTriangle(v1, v2, v3);
	return triangle;
}

std::unordered_map<std::string, Colour> parseMTL(std::string filename){
	std::unordered_map<std::string, Colour> m;
	std::string line;
	std::ifstream fin;
	std::string colourName;
	fin.open(filename);
	while (fin)
	{
		getline(fin, line);
		std::vector<std::string> tokenisedLine = split(line, ' ');
		if (tokenisedLine[0] == "newmtl") colourName = tokenisedLine[1];
		if (tokenisedLine[0] == "Kd") 
		{
			m[colourName] = Colour(round(std::stof(tokenisedLine[1])*255), round(std::stof(tokenisedLine[2])*255), round(std::stof(tokenisedLine[3])*255));
		}
	}
	fin.close();
	return m;
}

std::vector<ModelTriangle> parseOBJ(std::string fileName, float scale){
	std::vector<ModelTriangle> triangles;
	std::vector<glm::vec3> vertices;
	std::unordered_map<std::string, Colour> map;
	std::string line;
	std::string colourName;
	std::string mtlFileSource = "build/";
	std::ifstream fin;
	fin.open(fileName);
	while (fin)
	{
		getline(fin, line);
		std::vector<std::string> tokenisedLine = split(line, ' ');
		if(tokenisedLine[0] == "mtllib")
		{
			map = parseMTL(mtlFileSource.append(tokenisedLine[1]));
		}
		if (tokenisedLine[0] == "usemtl") colourName = tokenisedLine[1];
		if (tokenisedLine[0] == "v")
		{
			glm::vec3 vertex = glm::vec3(scale*std::stof(tokenisedLine[1]), scale*std::stof(tokenisedLine[2]), scale*std::stof(tokenisedLine[3]));
			vertices.push_back(vertex);
		}
		if (tokenisedLine[0] == "vn")
		{

		}
		if (tokenisedLine[0] == "f")
		{
			std::vector<std::vector<std::string>> tokenisedLine2;
			for(int i = 1; i < tokenisedLine.size(); i++){
				std::vector<std::string> tokenisedToken = split(tokenisedLine[i], '/');
				tokenisedLine2.push_back(tokenisedToken);
			}
			int v1 = std::stoi(tokenisedLine2[0][0]);
			int v2 = std::stoi(tokenisedLine2[1][0]);
			int v3 = std::stoi(tokenisedLine2[2][0]);
			ModelTriangle triangle = ModelTriangle(vertices[(v1-1)], vertices[(v2-1)], vertices[(v3-1)], map[colourName]);
			glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
			glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
			triangle.normal = glm::normalize(glm::cross(e0, e1));
			triangles.push_back(triangle);

		}
	}
	fin.close();
	return triangles;
	
}

CanvasPoint getCanvasIntersectionPoint( glm::vec3 vertexPosition, float focalLength, float depthArray[WIDTH][HEIGHT]) {
	CanvasPoint vertex;
	glm::vec3 vPosTransposed = vertexPosition - cameraPosition;
	vPosTransposed = vPosTransposed*cameraOrientation;
	float ratio = (focalLength)/vPosTransposed[2];
	vertex.x = round(-((ratio * vPosTransposed[0])*600) + (WIDTH/2));
	vertex.y = round(((ratio * vPosTransposed[1])*600) + (HEIGHT/2));
	vertex.depth = vPosTransposed.z;
	return vertex;
}

void drawObject(DrawingWindow &window, CanvasTriangle triangle, Colour colourGiven, float depthArray[WIDTH][HEIGHT]) {
	uint32_t colour = (255 << 24) + (colourGiven.red << 16) + (colourGiven.green << 8) + colourGiven.blue;
	std::vector<CanvasPoint> sortedVertices = sortVertices(triangle);
	CanvasPoint midPoint = calculateMidpoint(sortedVertices);  

	std::vector<CanvasPoint> topLine1 = interpolateLine(sortedVertices[0], midPoint);
	std::vector<CanvasPoint> topLine2 = interpolateLine(sortedVertices[0], sortedVertices[1]);
	float size1 = topLine1.size();
	float size2 = topLine2.size();
	//call the rake drawing function, always from left to right.
	if(size1 > 0 && size2 > 0){
		if(topLine1[round(size1/2)].x <= topLine2[round(size2/2)].x){
			for(uint32_t i = 0; i < size1; i++) {
				for(uint32_t j = 0; j < size2; j++) {
					if(topLine1[i].y == topLine2[j].y){
						drawObjectLine(window, topLine2[j], topLine1[i], colour, depthArray);
					}
				}
			}		
		}else{
			for(uint32_t i = 0; i < size1; i++) {
				for(uint32_t j = 0; j < size2; j++) {
					if(topLine1[i].y == topLine2[j].y){
						drawObjectLine(window, topLine1[i], topLine2[j], colour, depthArray);
					}
				}
			}
		}
	}
	
	std::vector<CanvasPoint> bottomLine1 = interpolateLine(midPoint, sortedVertices[2]);
	std::vector<CanvasPoint> bottomLine2 = interpolateLine(sortedVertices[1], sortedVertices[2]);
	if(bottomLine1.size() > 0 && bottomLine2.size() > 0){
		if(bottomLine1[bottomLine1.size()/2].x < bottomLine2[bottomLine2.size()/2].x){
			for(uint32_t i = 0; i < bottomLine1.size(); i++) {
				for(uint32_t j = 0; j < bottomLine2.size(); j++) {
					if (bottomLine2[j].y == bottomLine1[i].y){
						drawObjectLine(window, bottomLine2[j], bottomLine1[i], colour, depthArray);
					}
				}	
			}
		}else{
			for(uint32_t i = 0; i < bottomLine1.size(); i++) {
				for(uint32_t j = 0; j < bottomLine2.size(); j++) {
					if (bottomLine2[j].y == bottomLine1[i].y){
						drawObjectLine(window, bottomLine1[i], bottomLine2[j], colour, depthArray);
					}
				}
			}		
		}
	}
}

void lookAt(glm::vec3 focusPoint) {
	glm::vec3 forward = glm::normalize(cameraPosition - focusPoint);
	glm::vec3 right = glm::normalize(-glm::cross(forward, glm::vec3(0,1,0)));
	glm::vec3 up = glm::normalize(-glm::cross(right, forward));
	cameraOrientation = glm::mat3(right, up, forward);
}

RayTriangleIntersection  getClosestIntersection(glm::vec3 rayDirection, std::vector<ModelTriangle> triangles, glm::vec3 raySource){
	RayTriangleIntersection intersect;
	intersect.distanceFromCamera = 0;
	glm::vec3 edge1;
	glm::vec3 edge2;
 	for(float i=0; i < triangles.size(); i++){
		ModelTriangle triangle = triangles[i];
		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = raySource - triangle.vertices[0];
		glm::mat3 DEMatrix(-rayDirection, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
		if(possibleSolution[1] >= 0 && possibleSolution[1] <= 1) {
			if(possibleSolution[2] >= 0 && possibleSolution[2] <= 1) {
				if((possibleSolution[1] + possibleSolution[2]) <= 1) {
					//if no intersection has been found yet:
					if(intersect.distanceFromCamera == 0){
						intersect.distanceFromCamera = possibleSolution[0];
						intersect.intersectedTriangle = triangle;
						intersect.triangleIndex = i;
						//calculates the coordinates of the point with respect to the origin of the scence.
						glm::vec3 point= triangle.vertices[0] + (e0*possibleSolution[1]) + (e1*possibleSolution[2]);
						intersect.intersectionPoint = point;
					}
					//if an intersection has already been found -> check if closer than other intersection and replace if so.
					else if(possibleSolution[0] < intersect.distanceFromCamera && possibleSolution[0] > 0) {
						intersect.distanceFromCamera = possibleSolution[0];
						intersect.intersectedTriangle = triangle;
						intersect.triangleIndex = i;
						glm::vec3 point= triangle.vertices[0] + (e0*possibleSolution[1]) + (e1*possibleSolution[2]);
						intersect.intersectionPoint = point;
					}
				}
			}
		}
	}
	return intersect;
}

RayTriangleIntersection  getClosestShadowIntersection(glm::vec3 rayDirection, std::vector<ModelTriangle> triangles, glm::vec3 raySource, float triangleIndex){
	RayTriangleIntersection intersect;
	intersect.distanceFromCamera = 0;
 	for(float i=0; i < triangles.size(); i++){
		if(i != triangleIndex){
			ModelTriangle triangle = triangles[i];
			glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
			glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
			glm::vec3 SPVector = raySource - triangle.vertices[0];
			glm::mat3 DEMatrix(-rayDirection, e0, e1);
			glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;
			if(possibleSolution[1] >= 0 && possibleSolution[1] <= 1) {
				if(possibleSolution[2] >= 0 && possibleSolution[2] <= 1) {
					if((possibleSolution[1] + possibleSolution[2]) <= 1) {
						//if no intersection has been found yet:
						if(intersect.distanceFromCamera == 0){
							intersect.distanceFromCamera = possibleSolution[0];
							intersect.intersectedTriangle = triangle;
							intersect.triangleIndex = i;
							//calculates the coordinates of the point with respect to the origin of the scence.
							glm::vec3 point= triangle.vertices[0] + (e0*possibleSolution[1]) + (e1*possibleSolution[2]);
							intersect.intersectionPoint = point;
						}
						//if an intersection has already been found -> check if closer than other intersection and replace if so.
						else if(possibleSolution[0] < intersect.distanceFromCamera && possibleSolution[0] > 0) {
							intersect.distanceFromCamera = possibleSolution[0];
							intersect.intersectedTriangle = triangle;
							intersect.triangleIndex = i;
							glm::vec3 point= triangle.vertices[0] + (e0*possibleSolution[1]) + (e1*possibleSolution[2]);
							intersect.intersectionPoint = point;
						}
					}
				}
			}
		}
	}
	return intersect;
}

bool isPointInShadow(glm::vec3 point, std::vector<ModelTriangle> triangles, float triangleIndex) {
	glm::vec3 shadowRay = lightSource - point;
	RayTriangleIntersection intersect = getClosestShadowIntersection(glm::normalize(shadowRay), triangles, point, triangleIndex);
	//if no intersection found
	if(intersect.distanceFromCamera <= 0) return false;
	float distShadowRay = sqrt(pow(shadowRay.x,2) + pow(shadowRay.y,2) + pow(shadowRay.z,2));
	glm::vec3 intersectVec = intersect.intersectionPoint - point;
	float distIntersect = sqrt(pow(intersectVec.x,2) + pow(intersectVec.y,2) + pow(intersectVec.z,2));
	if(distIntersect < distShadowRay){
		return true;
	}
	else {
		return false;
	}
}

//function calculate the brightness scalar for specular lighting
float specular(ModelTriangle triangle, RayTriangleIntersection intersect){
	float brightness;
	glm::vec3 intersectToCamera = glm::normalize(cameraPosition - intersect.intersectionPoint);
	glm::vec3 lightToIntersect = glm::normalize(intersect.intersectionPoint - lightSource);
	glm::vec3 reflectedRay = glm::normalize(lightToIntersect - (2*triangle.normal*glm::dot(lightToIntersect, triangle.normal)));
	brightness = pow(glm::dot(intersectToCamera, reflectedRay), 8);
	if(brightness < 0) brightness = 0;
	return brightness;
}

float angleOfIncidence(ModelTriangle triangle, RayTriangleIntersection intersect, glm::vec3 intersectToLight){
	glm::vec3 intersectToLightNormalised  = glm::normalize(intersectToLight);
	float angleOfIncidence = glm::dot(triangle.normal, intersectToLightNormalised);
	float AoIBrightness = angleOfIncidence;
	if(AoIBrightness < 0) AoIBrightness = 0;
	return AoIBrightness;

}

float proximity(RayTriangleIntersection intersect ){
	float distToLight = glm::length(lightSource - intersect.intersectionPoint);
	float proxBrightness = 1/(6*3.1415*pow(distToLight, 2));
	if (proxBrightness > 1) proxBrightness = 1;
	return proxBrightness;
}


void draw(DrawingWindow &window, std::vector<ModelTriangle> triangles, float focalLength) {
	window.clearPixels();
	if (orbit == true) {
		glm::mat3 rotationMatrix = glm::mat3(glm::vec3(cosf(0.02), 0, -sinf(0.02)), glm::vec3(0, 1, 0), glm::vec3(sinf(0.02), 0, cosf(0.02)));
		cameraPosition = cameraPosition*rotationMatrix;
		lookAt(glm::vec3(0,0,0));
	}	
	for(float y = 0; y < HEIGHT; y++){
		for(float x = 0; x < WIDTH; x++){
			glm::vec3 rayDirection;
			glm::vec3 cameraSpacePixel;
			cameraSpacePixel.x = (x - WIDTH/2)/(WIDTH);
			cameraSpacePixel.y = (HEIGHT/2 - y)/(WIDTH);
			cameraSpacePixel.z = -focalLength;
 			rayDirection = glm::normalize(cameraSpacePixel*glm::inverse(cameraOrientation));
			
			RayTriangleIntersection intersect = getClosestIntersection(rayDirection, triangles, cameraPosition);
			if(intersect.distanceFromCamera > 0) {
				ModelTriangle triangle = triangles[intersect.triangleIndex];
				if(isPointInShadow(intersect.intersectionPoint, triangles, intersect.triangleIndex) == false){
					glm::vec3 intersectToLight = lightSource - intersect.intersectionPoint;
					float proxBrightness = proximity(intersect);
					float AoIBrightness = angleOfIncidence(triangle, intersect, intersectToLight);
					float specularBrightness = specular(triangle, intersect);
					float brightness = (AoIBrightness + proxBrightness + specularBrightness)/3;
					brightness = brightness + 0.2;
					if(brightness > 1) brightness =1;
					uint32_t colour = (255 << 24) + (static_cast<int>(round(triangle.colour.red * brightness)) << 16) + (static_cast<int>(round(triangle.colour.green*brightness)) << 8) + static_cast<int>(round(triangle.colour.blue*brightness));
					window.setPixelColour(x, y, colour);
				}
				else {
					float brightness = 0.2;
					uint32_t colour = (255 << 24) + (static_cast<int>(round(triangle.colour.red * brightness)) << 16) + (static_cast<int>(round(triangle.colour.green*brightness)) << 8) + static_cast<int>(round(triangle.colour.blue*brightness));
					window.setPixelColour(x, y, colour);
				}
			}
		}
	}
}

void drawRasterised(DrawingWindow &window, std::vector<ModelTriangle> triangles, float focalLength){
	window.clearPixels();
	if (orbit == true) {
		glm::mat3 rotationMatrix = glm::mat3(glm::vec3(cosf(0.02), 0, -sinf(0.02)), glm::vec3(0, 1, 0), glm::vec3(sinf(0.02), 0, cosf(0.02)));
		cameraPosition = cameraPosition*rotationMatrix;
		lookAt(glm::vec3(0,0,0));
	}
	float depthArray[WIDTH][HEIGHT] = {{0}};
	for(float i = 0; i < triangles.size(); i++){
		CanvasPoint v0 = getCanvasIntersectionPoint(triangles[i].vertices[0], focalLength, depthArray);
		CanvasPoint v1 = getCanvasIntersectionPoint( triangles[i].vertices[1], focalLength, depthArray);
		CanvasPoint v2 = getCanvasIntersectionPoint( triangles[i].vertices[2], focalLength, depthArray);
		CanvasTriangle triangle = CanvasTriangle(v0,v1,v2); 
		drawObject(window, triangle, triangles[i].colour, depthArray);
	}
}

void drawWireFrame(DrawingWindow &window, std::vector<ModelTriangle> triangles, float focalLength){
	window.clearPixels();
	if (orbit == true) {
		glm::mat3 rotationMatrix = glm::mat3(glm::vec3(cosf(0.02), 0, -sinf(0.02)), glm::vec3(0, 1, 0), glm::vec3(sinf(0.02), 0, cosf(0.02)));
		cameraPosition = cameraPosition*rotationMatrix;
		lookAt(glm::vec3(0,0,0));
	}
	float depthArray[WIDTH][HEIGHT] = {{0}};
	for(float i = 0; i < triangles.size(); i++){
		CanvasPoint v0 = getCanvasIntersectionPoint(triangles[i].vertices[0], focalLength, depthArray);
		CanvasPoint v1 = getCanvasIntersectionPoint( triangles[i].vertices[1], focalLength, depthArray);
		CanvasPoint v2 = getCanvasIntersectionPoint( triangles[i].vertices[2], focalLength, depthArray);
		CanvasTriangle triangle = CanvasTriangle(v0,v1,v2); 
		drawTriangle(window, triangle, triangles[i].colour);
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) cameraPosition = cameraPosition - glm::vec3(0.03, 0, 0);
		else if (event.key.keysym.sym == SDLK_RIGHT) cameraPosition = cameraPosition + glm::vec3(0.03, 0, 0);
		else if (event.key.keysym.sym == SDLK_DOWN) cameraPosition = cameraPosition - glm::vec3(0, 0.02, 0);
		else if (event.key.keysym.sym == SDLK_UP) cameraPosition = cameraPosition + glm::vec3(0, 0.02, 0);
		else if (event.key.keysym.sym == SDLK_z) cameraPosition = cameraPosition - glm::vec3(0,0, 0.1);
		else if (event.key.keysym.sym == SDLK_x) cameraPosition = cameraPosition + glm::vec3(0,0, 0.1);
		else if (event.key.keysym.sym == SDLK_o) {
			if(orbit == true) {
				orbit = false;
				std::cout << "stopping orbit" << std::endl;
			}
			else {
				std::cout << "starting orbit" << std::endl;
				orbit = true;
			}
		}
		else if (event.key.keysym.sym == SDLK_j) {

			glm::mat3 rotationMatrix = glm::mat3(glm::vec3(cosf(0.0174533), 0, -sinf(0.0174533)), glm::vec3(0, 1, 0), glm::vec3(sinf(0.0174533), 0, cosf(0.0174533)));
			cameraPosition = cameraPosition*rotationMatrix;
			lookAt(glm::vec3(0,0,0));
		}
		else if (event.key.keysym.sym == SDLK_l) {
			glm::vec3 col1 = glm::vec3(cosf(-0.0174533), 0, -sinf(-0.0174533));
			glm::vec3 col2 = glm::vec3(0, 1, 0);
			glm::vec3 col3 = glm::vec3(sinf(-0.0174533), 0, cosf(-0.0174533));
			glm::mat3 rotationMatrix = glm::mat3(col1, col2, col3);
			cameraPosition = cameraPosition*rotationMatrix;
			lookAt(glm::vec3(0,0,0));
		}
		else if(event.key.keysym.sym == SDLK_k){
			glm::vec3 col1 = glm::vec3(1, 0, 0);
			glm::vec3 col2 = glm::vec3(0, cosf(0.0174533), sinf(0.0174533));
			glm::vec3 col3 = glm::vec3(0, -sinf(0.0174533), cosf(0.0174533));
			glm::mat3 rotationMatrix = glm::mat3(col1, col2, col3);
			cameraPosition = cameraPosition*rotationMatrix;
			lookAt(glm::vec3(0,0,0));
		}
		else if(event.key.keysym.sym == SDLK_i){
			glm::vec3 col1 = glm::vec3(1, 0, 0);
			glm::vec3 col2 = glm::vec3(0, cosf(-0.0174533), sinf(-0.0174533));
			glm::vec3 col3 = glm::vec3(0, -sinf(-0.0174533), cosf(-0.0174533));
			glm::mat3 rotationMatrix = glm::mat3(col1, col2, col3);
			cameraPosition = cameraPosition*rotationMatrix;
			lookAt(glm::vec3(0,0,0));
		}
		else if (event.key.keysym.sym == SDLK_a) {
			glm::vec3 col1 = glm::vec3(cosf(-0.0174533), 0, -sinf(-0.0174533));
			glm::vec3 col2 = glm::vec3(0, 1, 0);
			glm::vec3 col3 = glm::vec3(sinf(-0.0174533), 0, cosf(-0.0174533));
			glm::mat3 rotationMatrix = glm::mat3(col1, col2, col3);
			cameraOrientation = cameraOrientation*rotationMatrix;
		}
		else if (event.key.keysym.sym == SDLK_d) {
			glm::mat3 rotationMatrix = glm::mat3(glm::vec3(cosf(0.0174533), 0, -sinf(0.0174533)), glm::vec3(0, 1, 0), glm::vec3(sinf(0.0174533), 0, cosf(0.0174533)));
			cameraOrientation = cameraOrientation*rotationMatrix;
		}
		else if(event.key.keysym.sym == SDLK_w){
			glm::mat3 rotationMatrix = glm::mat3(glm::vec3(1, 0, 0), glm::vec3(0, cosf(0.0174533), sinf(0.0174533)), glm::vec3(0, -sinf(0.0174533), cosf(0.0174533)));
			cameraOrientation = cameraOrientation*rotationMatrix;
		}
		else if(event.key.keysym.sym == SDLK_s){
			glm::mat3 rotationMatrix = glm::mat3(glm::vec3(1, 0, 0), glm::vec3(0, cosf(-0.0174533), sinf(-0.0174533)), glm::vec3(0, -sinf(-0.0174533), cosf(-0.0174533)));
			cameraOrientation = cameraOrientation*rotationMatrix;
		}			
		else if (event.key.keysym.sym == SDLK_u) {
			CanvasTriangle triangle = randTriangle(); 
			Colour colour = Colour(rand()%256,rand()%256,rand()%256);
			drawTriangle(window, triangle, colour);
		}
		else if (event.key.keysym.sym == SDLK_f) {
		
			CanvasTriangle triangle = randTriangle();
			Colour colour = Colour(rand()%256,rand()%256,rand()%256);
			drawFilledTriangle(window, triangle, colour);
		}
		else if (event.key.keysym.sym == SDLK_1) {
			std::cout << "Ray Traced Rendering" << std::endl;
			renderType = 1;
		}
		else if (event.key.keysym.sym == SDLK_2) {
			std::cout << "Rasterised Rendering" << std::endl;
			renderType = 2;
		}
		else if (event.key.keysym.sym == SDLK_3) {
			std::cout << "WireFrame Rendering" << std::endl;
			renderType = 3;
		}
		else if (event.key.keysym.sym == SDLK_MINUS) {
			lightSource.z = lightSource.z + 0.02;
		}
		else if (event.key.keysym.sym == SDLK_EQUALS) {
			lightSource.z = lightSource.z - 0.02;
		}	
		else if (event.key.keysym.sym == SDLK_PAGEDOWN) {
			lightSource.y = lightSource.y - 0.02;
		}	
		else if (event.key.keysym.sym == SDLK_PAGEUP) {
			lightSource.y = lightSource.y + 0.02;
		}
		else if (event.key.keysym.sym == SDLK_0) {
			lightSource.x = lightSource.x + 0.02;
		}	
		else if (event.key.keysym.sym == SDLK_9) {
			lightSource.x = lightSource.x - 0.02;
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	std::vector<ModelTriangle> triangles;	
	float scale = 0.17;
	triangles = parseOBJ("build/cornell-box.obj", scale);
	float focalLength = 2.0;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);		
		if(renderType == 1) draw(window, triangles, focalLength);
		else if(renderType == 2) drawRasterised(window, triangles, focalLength);
		else if(renderType == 3) drawWireFrame(window, triangles, focalLength);
		// std::cout << "(" << cameraPosition.x << "," << cameraPosition.y << ", " << cameraPosition.z << ")"  << std::endl;
		// Need to render the frame at the end, or nothing actually gets shown on the screen !		
		window.renderFrame();
	}
}