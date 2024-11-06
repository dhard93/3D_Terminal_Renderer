#include <iostream>
#include <sys/ioctl.h>
#include <unistd.h>
#include <vector>
#include <cmath>
#include <chrono>
#include <thread>
#include <algorithm>

// Todo:
// ! Figure out why line drawing not working right
// ! Add rotate method that can rotate along x, y, or z axis
// ! Add timer and controls to start/stop rotation
// ! Figure out how to scale x values to make up for fact that there is more vertical space in terminal
// ! Change Triangle and Mesh to classes
// ! Add class for matrix/vector operations
// ! Add fill method to fill triangles using ASCII chars and scanlines
// ! Implement z-buffer to only render what is visible to camera
// ! Implement light source that determines what chars are rendered based on grayscale values
// ! Add different 3d shapes that mesh class can be used to represent
// ! Add controls to add/remove shapes and apply basic physics to them.
// ! Add control to roll like dice and 

// After tutorial, final test is to implement this project as 'Demo' mode in imgProcessor web application.
// It should allow user to draw 3d shapes, rotate/move them, add/remove shapes from canvas, apply physics to them,
// fill them with various colors/characters using different techniques such as scanlines, solid fill, ASCII art, etc...
// Inspired by demoscene where canvas is split in half and each half places different visual effects on the shapes,
// depending which side they are on.
// Use to create dice roller that has dice textures mapped to each cube face.

// Write down all math involved for projection, rotation, etc... and then implement on own without any help.
// Go through the math and work it out until you have an actual understanding of why it is the way it is.

struct vec3
{
    float x, y, z;
};

struct vec4
{
    float w, x, y, z;
};

struct edge
{
    vec3 start;
    vec3 end;
};

struct triangle
{
    vec3 vertices[3];
    edge edges[3];

    triangle() {};

    triangle(vec3 a, vec3 b, vec3 c)
    {
        vertices[0] = a;
        vertices[1] = b;
        vertices[2] = c;
        edges[0] = {a, b};
        edges[1] = {b, c};
        edges[2] = {a, c};
    }
};

struct mesh
{
    std::vector<triangle> triangles;
};

struct mat3x3
{
    vec3 rows[3];
};

struct mat4x4
{
    vec4 rows[4];
};

void getTerminalSize(int &width, int &height);
void multiplyVecMatrix(vec3 &inPoint, vec3 &outPoint, const mat4x4 &mat);
void drawPoint(const float &xValue, const float &yValue, const int &width, const int &height) ;
void drawTriangles(const mesh &cubeMesh, const int width, const int height, const float theta, const mat4x4 &projMat);
void rotate(const char axis, const vec3 &inVec, vec3 &outVec, const float &theta);
void fillTriangle(triangle &tri, const float &width, const float &height, const float &dir);
float interpolate(const float minA, const float minB, const float maxA, const float maxB, const int cur);
void sort(vec3 verts[3], const int code);
void drawChar(const int &sWidth, const int &sHeight, const int &xVal, const int &yVal, const float &zVal);

int main()
{
    int screenW, screenH;
    float aspectRatio;
    float fov = 90.0f;
    float fNear = 0.1f;
    float fFar = 1000.0f;
    float fovRad = 1.0f/tanf(fov * 0.5f / 180.0f * 3.14159f);

    std::cout << "\033[?25l"; // Hide the cursor
    getTerminalSize(screenW, screenH);
    aspectRatio = (float)screenH/(float)screenW;

    vec3 pointA = {0, 0, 0};
    vec3 pointB = {0, 0, 1};
    vec3 pointC = {0, 1, 0};
    vec3 pointD = {0, 1, 1};
    vec3 pointE = {1, 0, 0};
    vec3 pointF = {1, 0, 1};
    vec3 pointG = {1, 1, 0};
    vec3 pointH = {1, 1, 1};

    mesh cubeMesh;

    cubeMesh.triangles =
    {
        // South Face
        triangle(pointA, pointC, pointG),
        triangle(pointA, pointG, pointE),

        // East Face
        triangle(pointE, pointG, pointH),
        triangle(pointE, pointH, pointF),

        // North Face
        triangle(pointF, pointH, pointD),
        triangle(pointF, pointD, pointB),

        // West Face
        triangle(pointB, pointD, pointC),
        triangle(pointB, pointC, pointA),

        // Top Face
        triangle(pointC, pointD, pointH),
        triangle(pointC, pointH, pointG),

        // Bottom Face
        triangle(pointF, pointB, pointA),
        triangle(pointF, pointA, pointE)
    };

    float fTheta = fovRad;

    // Projection Matrix
    mat4x4 projMat = {};
    projMat.rows[0].w = aspectRatio * fovRad;
    projMat.rows[1].x = fovRad;
    projMat.rows[2].y = (fFar / (fFar - fNear));
    projMat.rows[3].y = ((-fFar * fNear) / (fFar - fNear));
    projMat.rows[2].z = 1.0f;
    projMat.rows[3].z = 0.0f;

    mat4x4 rotMatX = {};
    rotMatX.rows[0].w = 1.0f;
    rotMatX.rows[1].x = cosf(fTheta);
    rotMatX.rows[1].y = -sinf(fTheta);
    rotMatX.rows[2].x = sinf(fTheta);
    rotMatX.rows[2].y = cosf(fTheta);
    rotMatX.rows[3].z = 1.0f;

    mat4x4 rotMatY = {};
    rotMatY.rows[0].w = cosf(fTheta * 0.5f);
    rotMatY.rows[0].y = sinf(fTheta * 0.5f);
    rotMatY.rows[1].y = 1.0f;
    rotMatY.rows[2].w = -sinf(fTheta * 0.5f);
    rotMatY.rows[2].y = cosf(fTheta * 0.5f);
    rotMatY.rows[3].z = 1.0f;

    mat4x4 rotMatZ = {};
    rotMatZ.rows[0].w = cosf(fTheta * 0.25f);
    rotMatZ.rows[0].x = -sinf(fTheta * 0.25f);
    rotMatZ.rows[1].w = sinf(fTheta * 0.25f);
    rotMatZ.rows[1].x = cosf(fTheta * 0.25f);
    rotMatZ.rows[2].y = 1.0f;
    rotMatZ.rows[3].z = 1.0f;

    auto lastTime = std::chrono::high_resolution_clock::now();

    while (true)
    {
        auto currentTime = std::chrono::high_resolution_clock::now();
        auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - lastTime);

        if (elapsedTime.count() >= 100)
        {
            system("clear");
            drawTriangles(cubeMesh, screenW, screenH, fTheta, projMat);
            std::cout << std::flush;
            fTheta += .1;
            lastTime = currentTime;
        }
    }

    // On quit
    std::cout << "\033[?25h"; // Show the cursor
    system("clear");
    std::cout << "\033[" << screenH << ";1H";

    return 0;
}

void getTerminalSize(int &width, int &height)
{
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    width = w.ws_col;
    height = w.ws_row;
    system("clear");
}

// Multply Point by projection matrix to get projected Point
void multiplyVecMatrix(vec3 &inPoint, vec3 &outPoint, const mat4x4 &mat)
{
    float w = inPoint.x * mat.rows[0].z + inPoint.y * mat.rows[1].z + inPoint.z * mat.rows[2].z + mat.rows[3].z;

    if (w != 0)
    {
        outPoint.x = (inPoint.x * mat.rows[0].w * 2 + inPoint.y * mat.rows[1].w + inPoint.z * mat.rows[2].w + mat.rows[3].w) / w;
        outPoint.y = (inPoint.x * mat.rows[0].x * 2 + inPoint.y * mat.rows[1].x + inPoint.z * mat.rows[2].x + mat.rows[3].x) / w;
        outPoint.z = (inPoint.x * mat.rows[0].y * 2 + inPoint.y * mat.rows[1].y + inPoint.z * mat.rows[2].y + mat.rows[3].y) / w;
    }
}

void drawPoint(const float &xValue, const float &yValue, const int &width, const int &height) {
    int x = floor(xValue);
    int y = floor(yValue);

    // Ensure the point is within terminal bounds
    if (x >= 0 && x <= width && y >= 0 && y <= height) {
        // Move the cursor to the specified point
        std::cout << "\033[" << y << ";" << x << "H*";
    }
}

void rotate(const char axis, const vec3 &inVec, vec3 &outVec, const float &theta)
{
    mat3x3 rotMat {};

    switch (axis)
    {
        case 'x':
            rotMat.rows[0].x = 1;
            rotMat.rows[1].y = cos(theta);
            rotMat.rows[1].z = -sin(theta);
            rotMat.rows[2].y = sin(theta);
            rotMat.rows[2].z = cos(theta);
            break;
        case 'y':
            rotMat.rows[0].x = cosf(theta);
            rotMat.rows[0].z = sinf(theta);
            rotMat.rows[1].y = 1.0f;
            rotMat.rows[2].x = -sinf(theta);
            rotMat.rows[2].z = cosf(theta);
            break;
        case 'z':
            rotMat.rows[0].x = cos(theta);
            rotMat.rows[0].y = -sin(theta);
            rotMat.rows[1].x = sin(theta);
            rotMat.rows[1].y = cos(theta);
            rotMat.rows[2].z = 1;
    }

    outVec.x = inVec.x * rotMat.rows[0].x + inVec.y * rotMat.rows[1].x + inVec.z * rotMat.rows[2].x;
    outVec.y = inVec.x * rotMat.rows[0].y + inVec.y * rotMat.rows[1].y + inVec.z * rotMat.rows[2].y;
    outVec.z = inVec.x * rotMat.rows[0].z + inVec.y * rotMat.rows[1].z + inVec.z * rotMat.rows[2].z;
}

void drawTriangles(const mesh &cubeMesh, const int width, const int height, const float theta, const mat4x4 &projMat)
{
    // Draw Triangles by projecting Each of their Points into screen space.
    for (auto tri : cubeMesh.triangles)
    {

        triangle translatedTri, projTri, rotatedTriX, rotatedTriZ;

        rotate('z', tri.vertices[0], rotatedTriZ.vertices[0], theta);
        rotate('z', tri.vertices[1], rotatedTriZ.vertices[1], theta);
        rotate('z', tri.vertices[2], rotatedTriZ.vertices[2], theta);

        rotate('x', rotatedTriZ.vertices[0], rotatedTriX.vertices[0], theta * 0.5f);
        rotate('x', rotatedTriZ.vertices[1], rotatedTriX.vertices[1], theta * 0.5f);
        rotate('x', rotatedTriZ.vertices[2], rotatedTriX.vertices[2], theta * 0.5f);

        translatedTri = rotatedTriX;
        translatedTri.vertices[0].z += 3.0f;
        translatedTri.vertices[1].z += 3.0f;
        translatedTri.vertices[2].z += 3.0f;

        vec3 AB = {translatedTri.vertices[1].x - translatedTri.vertices[0].x, translatedTri.vertices[1].y - translatedTri.vertices[0].y, translatedTri.vertices[1].z - translatedTri.vertices[0].z};
        vec3 AC = {translatedTri.vertices[2].x - translatedTri.vertices[0].x, translatedTri.vertices[2].y - translatedTri.vertices[0].y, translatedTri.vertices[2].z - translatedTri.vertices[0].z};
        vec3 normal = {(AB.y * AC.z) - (AB.z * AC.y), (AB.z * AC.x) - (AB.x * AC.z), (AB.x * AC.y) - (AB.y * AC.x)};
        float mag = sqrtf(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
        normal.x /= mag;
        normal.y /= mag;
        normal.z /= mag;

        if (normal.x * translatedTri.vertices[0].x + normal.y * translatedTri.vertices[0].y + normal.z * translatedTri.vertices[0].z < 0)
        {
            multiplyVecMatrix(translatedTri.vertices[0], projTri.vertices[0], projMat);
            multiplyVecMatrix(translatedTri.vertices[1], projTri.vertices[1], projMat);
            multiplyVecMatrix(translatedTri.vertices[2], projTri.vertices[2], projMat);

            projTri.vertices[0].x += 1.0f; projTri.vertices[0].y += 1.0f;
            projTri.vertices[1].x += 1.0f; projTri.vertices[1].y += 1.0f;
            projTri.vertices[2].x += 1.0f; projTri.vertices[2].y += 1.0f;

            projTri.vertices[0].x *= 0.5f * width;
            projTri.vertices[1].x *= 0.5f * width;
            projTri.vertices[2].x *= 0.5f * width;
            projTri.vertices[0].y *= 0.5f * height;
            projTri.vertices[1].y *= 0.5f * height;
            projTri.vertices[2].y *= 0.5f * height;

            //drawPoint(projTri.vertices[0].x, projTri.vertices[0].y, width, height);
            //drawPoint(projTri.vertices[1].x, projTri.vertices[1].y, width, height);
            //drawPoint(projTri.vertices[2].x, projTri.vertices[2].y, width, height);
            fillTriangle(projTri, width, height, normal.z);
        }
    }
}

// Fill triangle using scanline algorithm.
/*
    - sort triangle vertices by y value.
    - iterate triangle from m y value to max y value, 1 row at a time.
    - for each scanline, determine if it is above middle vertex or not:
        - if it is, then the scanline is iterating between collision points with edges AB and AC
        - if not, then the scanline is iterating between collision points with edges BC and AC.
    - using slopes of the 2 edges in use, interpolate x position of both edge collision points and sort in ascending order.
    - iterate from x0 to x1 of the scanline, interpolating z value of that position and drawing char based on z-depth.
*/
void fillTriangle(triangle &tri, const float &width, const float &height, const float &dir)
{
    // Construct bounding box and sort vertices from lowest to highest height.
    float minX = 0, maxX = 0;
    float minZ = 0, maxZ = 0;
    float z = 0;

    triangle sortedTri = tri;
    sort(sortedTri.vertices, 1);

    int minY = std::ceil(sortedTri.vertices[0].y);
    int maxY = std::floor(sortedTri.vertices[2].y);
    int temp;

    // Iterate from top to bottom y.
    // ! The issue is it is using y to interpolate min and max z when it should be using x.
    for (int y = minY; y <= maxY; ++y)
    {

        maxX = interpolate(sortedTri.vertices[0].x, sortedTri.vertices[0].y, sortedTri.vertices[2].x, sortedTri.vertices[2].y, y);

        if (y < sortedTri.vertices[1].y)
            minX =interpolate(sortedTri.vertices[0].x, sortedTri.vertices[0].y, sortedTri.vertices[1].x, sortedTri.vertices[1].y, y);
        else
            minX =interpolate(sortedTri.vertices[1].x, sortedTri.vertices[1].y, sortedTri.vertices[2].x, sortedTri.vertices[2].y, y);

        if (minX > maxX)
        {
            temp = minX;
            minX = maxX;
            maxX = temp;
        }

        // Interpolate z value of each char in scanline and choose char to
        // draw based on z depth.
        for (float x = minX; x <= maxX; ++x)
        {
            maxZ = interpolate(sortedTri.vertices[0].z, minX, sortedTri.vertices[2].z, maxX, x);
            if (y < sortedTri.vertices[1].y)
                minZ =interpolate(sortedTri.vertices[0].z, minX, sortedTri.vertices[1].z, maxX, x);
            else
                minZ =interpolate(sortedTri.vertices[1].z, minX, sortedTri.vertices[2].z, maxX, x);
            if (minZ > maxZ)
            {
                temp = minZ;
                minZ = maxZ;
                maxZ = temp;
            }
            //z = interpolate(minZ, minX, maxZ, maxX, x) * (-1) / dir;
            drawChar(width, height, x, y, dir);
            //std::cout << z << std::endl;
        }
    }
}

float interpolate(const float minA, const float minB, const float maxA, const float maxB, const int cur)
{ return minA + ((cur - minB) * (maxA - minA) / (maxB - minB)); }

void sort(vec3 verts[3], const int code)
{
    vec3 temp;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 1; j < 3; ++j)
        {
            if (code == 0)
            {
                if (verts[i].x > verts[j].x)
                {
                    temp = verts[j];
                    verts[j] = verts[i];
                    verts[i] = temp;
                }
            }
            else
            {
                if (verts[i].y > verts[j].y)
                {
                    temp = verts[j];
                    verts[j] = verts[i];
                    verts[i] = temp;
                }
            }
        }
    }
}

// ! Fix zVal thresholds !
void drawChar(const int &sWidth, const int &sHeight, const int &xVal, const int &yVal, const float &dir)
{
    char c;

    if (dir <=-.9)
        c = '#';
    else if (dir <=-0.8)
        c = '@';
    else if (dir <=-0.7)
        c = 'C';
    else if (dir <=-0.6)
        c = '>';
    else if (dir <=-0.5)
        c = 'c';
    else if (dir <=-0.4)
        c = '*';
    else if (dir <=-0.3)
        c = '-';
    else
        c = '.';

    // Ensure the point is within terminal bounds
    if (xVal >= 0 && xVal <= sWidth && yVal >= 0 && yVal <= sHeight) {
        // Move the cursor to the specified point
        std::cout << "\033[" << yVal << ";" << xVal << "H" << c;
    }
}