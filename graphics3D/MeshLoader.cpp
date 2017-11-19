
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(__APPLE__)
#include <GLUT/GLUT.h>
#include <OpenGL/gl3.h>
#include <OpenGL/glu.h>
#else
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#endif
#include <GL/glew.h>		 
#include <GL/freeglut.h>	
#endif

#include <string>
#include <vector>
#include <fstream>
#include <algorithm> 
const unsigned int windowWidth = 512, windowHeight = 512;

int majorVersion = 3, minorVersion = 0;

bool keyboardState[256];

void getErrorInfo(unsigned int handle) 
{
	int logLen;
	glGetShaderiv(handle, GL_INFO_LOG_LENGTH, &logLen);
	if (logLen > 0) 
	{
		char * log = new char[logLen];
		int written;
		glGetShaderInfoLog(handle, logLen, &written, log);
		printf("Shader log:\n%s", log);
		delete log;
	}
}

void checkShader(unsigned int shader, char * message) 
{
	int OK;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &OK);
	if (!OK) 
	{
		printf("%s!\n", message);
		getErrorInfo(shader);
	}
}

void checkLinking(unsigned int program) 
{
	int OK;
	glGetProgramiv(program, GL_LINK_STATUS, &OK);
	if (!OK) 
	{
		printf("Failed to link shader program!\n");
		getErrorInfo(program);
	}
}

// row-major matrix 4x4
struct mat4 
{
	float m[4][4];
public:
	mat4() {}
	mat4(float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33) 
	{
		m[0][0] = m00; m[0][1] = m01; m[0][2] = m02; m[0][3] = m03;
		m[1][0] = m10; m[1][1] = m11; m[1][2] = m12; m[1][3] = m13;
		m[2][0] = m20; m[2][1] = m21; m[2][2] = m22; m[2][3] = m23;
		m[3][0] = m30; m[3][1] = m31; m[3][2] = m32; m[3][3] = m33;
	}

	mat4 operator*(const mat4& right) 
	{
		mat4 result;
		for (int i = 0; i < 4; i++) 
		{
			for (int j = 0; j < 4; j++) 
			{
				result.m[i][j] = 0;
				for (int k = 0; k < 4; k++) result.m[i][j] += m[i][k] * right.m[k][j];
			}
		}
		return result;
	}
	operator float*() { return &m[0][0]; }
};


// 3D point in homogeneous coordinates
struct vec4 
{
	float v[4];

	vec4(float x = 0, float y = 0, float z = 0, float w = 1) 
	{
		v[0] = x; v[1] = y; v[2] = z; v[3] = w;
	}

	vec4 operator*(const mat4& mat) 
	{
		vec4 result;
		for (int j = 0; j < 4; j++) 
		{
			result.v[j] = 0;
			for (int i = 0; i < 4; i++) result.v[j] += v[i] * mat.m[i][j];
		}
		return result;
	}

	vec4 operator+(const vec4& vec) 
	{
		vec4 result(v[0] + vec.v[0], v[1] + vec.v[1], v[2] + vec.v[2], v[3] + vec.v[3]);
		return result;
	}
};

// 2D point in Cartesian coordinates
struct vec2 
{
	float x, y;

	vec2(float x = 0.0, float y = 0.0) : x(x), y(y) {}

	vec2 operator+(const vec2& v) 
	{
		return vec2(x + v.x, y + v.y);
	}

	vec2 operator*(float s) 
	{
		return vec2(x * s, y * s);
	}

};

// 3D point in Cartesian coordinates
struct vec3 
{
	float x, y, z;

	vec3(float x = 0.0, float y = 0.0, float z = 0.0) : x(x), y(y), z(z) {}

	static vec3 random() { return vec3(((float)rand() / RAND_MAX) * 2 - 1, ((float)rand() / RAND_MAX) * 2 - 1, ((float)rand() / RAND_MAX) * 2 - 1); }

	vec3 operator+(const vec3& v) { return vec3(x + v.x, y + v.y, z + v.z); }

	vec3 operator-(const vec3& v) { return vec3(x - v.x, y - v.y, z - v.z); }

	vec3 operator*(float s) { return vec3(x * s, y * s, z * s); }

	vec3 operator/(float s) { return vec3(x / s, y / s, z / s); }

	float length() { return sqrt(x * x + y * y + z * z); }

	vec3 normalize() { return *this / length(); }
    
	void print() { printf("%f \t %f \t %f \n", x, y, z); }
};

vec3 cross(const vec3& a, const vec3& b)
{
	return vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x );
}

class Wood{
    float scale;
    float turbulence;
    float period;
    float sharpness;
    
    Wood() 	{
        scale = 16;
        turbulence = 500;
        period = 8;
        sharpness = 10;
    }
    
    float f(vec3 r) {
        vec3 s = vec3(7502, 22777, 4767);
        float w = 0.0;
        for(int i=0; i<16; i++) {
            vec3 v1 = s - vec3(32768, 32768, 32768);
            vec3 v2 = r * 40.0;
            w += sin( (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z) / 65536.0);
            s = vec3((int)s.x % (int)32768.0 * 2.0,(int)s.y % (int)32768.0 * 2.0, (int)s.z % (int)32768.0 * 2.0) + vec3(floor(s.x / 32768.0),floor(s.y / 32768.0),floor(s.z / 32768.0));
        }
        return w / 32.0 + 0.5;
    }
    
    vec3 getColor(vec3 position) {
        float w =
        position.x * period +
        pow(f(position * scale), sharpness)*turbulence;
        w -= int(w + 10000.0);  // take fractional part
        return					// combine
        vec3(1, 0.3, 0) * w + 		// light wood
        vec3(0.35, 0.1, 0.05) * (1-w);	// dark wood
    }
};






class Geometry
{
protected:
	unsigned int vao;

public:
	Geometry()
	{
		glGenVertexArrays(1, &vao);					
	}

	virtual void Draw() = 0;
};


class TexturedQuad : public Geometry
{
    unsigned int vbo[3];
    
public:
    TexturedQuad()
    {
        glBindVertexArray(vao);
        
        glGenBuffers(3, &vbo[0]);
        
        //Vertex Position Array
        glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
        static float vertexCoords[] =
        {
            -1,  0, -1,
            -1,  0,  1,
             1,  0, -1,
             1,  0,  1 };
        
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertexCoords), vertexCoords, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);	// read 3 attributes at a time
        
        //Texture Coordinates Array
        glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
        static float textCoords[] =
        {
            0, 0,
            1, 0,
            0, 1,
            1, 1 };
        
        glBufferData(GL_ARRAY_BUFFER, sizeof(textCoords), textCoords, GL_STATIC_DRAW);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, NULL);
        
        
        //Normal Array
        glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
        static float normalCoords[] = {
            0, 1, 0,
            0, 1, 0,
            0, 1, 0,
            0, 1, 0
        };	// vertex data on the CPU
        
        glBufferData(GL_ARRAY_BUFFER,	// copy to the GPU
                     sizeof(normalCoords),	// size of the vbo in bytes
                     normalCoords,		// address of the data array on the CPU
                     GL_STATIC_DRAW);	// copy to that part of the memory which is not modified
        
        // map Attribute Array 0 to the currently bound vertex buffer (vbo)
        glEnableVertexAttribArray(2);
        
        // data organization of Attribute Array 2
        glVertexAttribPointer(2,	// Attribute Array 2
                              3, GL_FLOAT,		// components/attribute, component type
                              GL_FALSE,		// not in fixed point format, do not normalized
                              0, NULL);		// stride and offset: it is tightly packed
    }
    
    void Draw()
    {
        glEnable(GL_BLEND); // necessary for transparent pixels
        glEnable(GL_DEPTH_TEST);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glDisable(GL_BLEND);
        glDisable(GL_DEPTH_TEST);
    }
};

// inifnite ground

class InfiniteTexturedQuad : public Geometry {
    unsigned int vbo[3];
    
public:
    InfiniteTexturedQuad()
        {
            glBindVertexArray(vao);
            
            glGenBuffers(3, &vbo[0]);
            
            //Vertex Position Array
            glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
            static float vertexCoords[] =
            {
                0, 0, 0, 1,
                -1,  0, -1, 0,
                -1,  0,  1, 0,
                1,  0,  1,  0,
                1,  0,  -1, 0,
                -1, 0,  -1, 0
            };
            
            glBufferData(GL_ARRAY_BUFFER, sizeof(vertexCoords), vertexCoords, GL_STATIC_DRAW);
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, NULL);	// read 3 attributes at a time
            
            //Texture Coordinates Array
            glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
            static float textCoords[] =
            {
                0, 0,
                1, 0,
                0, 1,
                1, 1 };
            
            glBufferData(GL_ARRAY_BUFFER, sizeof(textCoords), textCoords, GL_STATIC_DRAW);
            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, NULL);
            
            
            //Normal Array
            glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
            static float normalCoords[] = {
                0, 1, 0,
                0, 1, 0,
                0, 1, 0,
                0, 1, 0
            };	// vertex data on the CPU
            
            glBufferData(GL_ARRAY_BUFFER,	// copy to the GPU
                         sizeof(normalCoords),	// size of the vbo in bytes
                         normalCoords,		// address of the data array on the CPU
                         GL_STATIC_DRAW);	// copy to that part of the memory which is not modified
            
            // map Attribute Array 0 to the currently bound vertex buffer (vbo)
            glEnableVertexAttribArray(2);
            
            // data organization of Attribute Array 2
            glVertexAttribPointer(2,	// Attribute Array 2
                                  3, GL_FLOAT,		// components/attribute, component type
                                  GL_FALSE,		// not in fixed point format, do not normalized
                                  0, NULL);		// stride and offset: it is tightly packed
        }
    
    void Draw() {
        glEnable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 6);
        glDisable(GL_BLEND);
        glDisable(GL_DEPTH_TEST);
    }
};

class  PolygonalMesh : public Geometry
{
	struct  Face
	{
		int       positionIndices[4];
		int       normalIndices[4];
		int       texcoordIndices[4];
		bool      isQuad;
	};

	std::vector<std::string*> rows;
	std::vector<vec3*> positions;
	std::vector<std::vector<Face*>> submeshFaces;
	std::vector<vec3*> normals;
	std::vector<vec2*> texcoords;

	int nTriangles;

public:
	PolygonalMesh(const char *filename);
	~PolygonalMesh();

	void Draw();
};



PolygonalMesh::PolygonalMesh(const char *filename)
{
	std::fstream file(filename); 
	if(!file.is_open())       
	{
		return;
	}

	char buffer[256];
	while(!file.eof())
	{
		file.getline(buffer,256);
		rows.push_back(new std::string(buffer));
	}

	submeshFaces.push_back(std::vector<Face*>());
	std::vector<Face*>* faces = &submeshFaces.at(submeshFaces.size()-1);

	for(int i = 0; i < rows.size(); i++)
	{
		if(rows[i]->empty() || (*rows[i])[0] == '#') 
			continue;      
		else if((*rows[i])[0] == 'v' && (*rows[i])[1] == ' ')
		{
			float tmpx,tmpy,tmpz;
			sscanf(rows[i]->c_str(), "v %f %f %f" ,&tmpx,&tmpy,&tmpz);      
			positions.push_back(new vec3(tmpx,tmpy,tmpz));  
		}
		else if((*rows[i])[0] == 'v' && (*rows[i])[1] == 'n')    
		{
			float tmpx,tmpy,tmpz;   
			sscanf(rows[i]->c_str(), "vn %f %f %f" ,&tmpx,&tmpy,&tmpz);
			normals.push_back(new vec3(tmpx,tmpy,tmpz));     
		}
		else if((*rows[i])[0] == 'v' && (*rows[i])[1] == 't')
		{
			float tmpx,tmpy;
			sscanf(rows[i]->c_str(), "vt %f %f" ,&tmpx,&tmpy);
			texcoords.push_back(new vec2(tmpx,tmpy));     
		}
		else if((*rows[i])[0] == 'f')  
		{
			if(count(rows[i]->begin(),rows[i]->end(), ' ') == 3)
			{
				Face* f = new Face();
				f->isQuad = false;
				sscanf(rows[i]->c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d",
					&f->positionIndices[0], &f->texcoordIndices[0], &f->normalIndices[0],
					&f->positionIndices[1], &f->texcoordIndices[1], &f->normalIndices[1],
					&f->positionIndices[2], &f->texcoordIndices[2], &f->normalIndices[2]);
				faces->push_back(f);
			}
			else
			{
				Face* f = new Face();
				f->isQuad = true;
				sscanf(rows[i]->c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d", 
					&f->positionIndices[0], &f->texcoordIndices[0], &f->normalIndices[0],
					&f->positionIndices[1], &f->texcoordIndices[1], &f->normalIndices[1],
					&f->positionIndices[2], &f->texcoordIndices[2], &f->normalIndices[2],
					&f->positionIndices[3], &f->texcoordIndices[3], &f->normalIndices[3]);
				faces->push_back(f);   
			}
		}
		else if((*rows[i])[0] == 'g')
		{
			if(faces->size() > 0)
			{
				submeshFaces.push_back(std::vector<Face*>());
				faces = &submeshFaces.at(submeshFaces.size()-1);
			}
		}
	}
	
	int numberOfTriangles = 0;
	for(int iSubmesh=0; iSubmesh<submeshFaces.size(); iSubmesh++)
	{
		std::vector<Face*>& faces = submeshFaces.at(iSubmesh);

		for(int i=0;i<faces.size();i++)
		{
			if(faces[i]->isQuad) numberOfTriangles += 2;
			else numberOfTriangles += 1;
		}
	}

	nTriangles = numberOfTriangles;
	
	float *vertexCoords = new float[numberOfTriangles * 9];
	float *vertexTexCoords = new float[numberOfTriangles * 6];
	float *vertexNormalCoords = new float[numberOfTriangles * 9];


	int triangleIndex = 0;
	for(int iSubmesh=0; iSubmesh<submeshFaces.size(); iSubmesh++)
	{
		std::vector<Face*>& faces = submeshFaces.at(iSubmesh);

		for(int i=0;i<faces.size();i++)
		{
			if(faces[i]->isQuad) 
			{
				vertexTexCoords[triangleIndex * 6] =     texcoords[faces[i]->texcoordIndices[0]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 1] = 1-texcoords[faces[i]->texcoordIndices[0]-1]->y;
				
				vertexTexCoords[triangleIndex * 6 + 2] = texcoords[faces[i]->texcoordIndices[1]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 3] = 1-texcoords[faces[i]->texcoordIndices[1]-1]->y;

				vertexTexCoords[triangleIndex * 6 + 4] = texcoords[faces[i]->texcoordIndices[2]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 5] = 1-texcoords[faces[i]->texcoordIndices[2]-1]->y;


				vertexCoords[triangleIndex * 9] =     positions[faces[i]->positionIndices[0]-1]->x;
				vertexCoords[triangleIndex * 9 + 1] = positions[faces[i]->positionIndices[0]-1]->y;
				vertexCoords[triangleIndex * 9 + 2] = positions[faces[i]->positionIndices[0]-1]->z;

				vertexCoords[triangleIndex * 9 + 3] = positions[faces[i]->positionIndices[1]-1]->x;
				vertexCoords[triangleIndex * 9 + 4] = positions[faces[i]->positionIndices[1]-1]->y;
				vertexCoords[triangleIndex * 9 + 5] = positions[faces[i]->positionIndices[1]-1]->z;

				vertexCoords[triangleIndex * 9 + 6] = positions[faces[i]->positionIndices[2]-1]->x;
				vertexCoords[triangleIndex * 9 + 7] = positions[faces[i]->positionIndices[2]-1]->y;
				vertexCoords[triangleIndex * 9 + 8] = positions[faces[i]->positionIndices[2]-1]->z;


				vertexNormalCoords[triangleIndex * 9] =     normals[faces[i]->normalIndices[0]-1]->x;    
				vertexNormalCoords[triangleIndex * 9 + 1] = normals[faces[i]->normalIndices[0]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 2] = normals[faces[i]->normalIndices[0]-1]->z;

				vertexNormalCoords[triangleIndex * 9 + 3] = normals[faces[i]->normalIndices[1]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 4] = normals[faces[i]->normalIndices[1]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 5] = normals[faces[i]->normalIndices[1]-1]->z;

				vertexNormalCoords[triangleIndex * 9 + 6] = normals[faces[i]->normalIndices[2]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 7] = normals[faces[i]->normalIndices[2]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 8] = normals[faces[i]->normalIndices[2]-1]->z;
								
				triangleIndex++;


				vertexTexCoords[triangleIndex * 6] =     texcoords[faces[i]->texcoordIndices[1]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 1] = 1-texcoords[faces[i]->texcoordIndices[1]-1]->y;
				
				vertexTexCoords[triangleIndex * 6 + 2] = texcoords[faces[i]->texcoordIndices[2]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 3] = 1-texcoords[faces[i]->texcoordIndices[2]-1]->y;

				vertexTexCoords[triangleIndex * 6 + 4] = texcoords[faces[i]->texcoordIndices[3]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 5] = 1-texcoords[faces[i]->texcoordIndices[3]-1]->y;


				vertexCoords[triangleIndex * 9] =     positions[faces[i]->positionIndices[1]-1]->x;
				vertexCoords[triangleIndex * 9 + 1] = positions[faces[i]->positionIndices[1]-1]->y;
				vertexCoords[triangleIndex * 9 + 2] = positions[faces[i]->positionIndices[1]-1]->z;

				vertexCoords[triangleIndex * 9 + 3] = positions[faces[i]->positionIndices[2]-1]->x;
				vertexCoords[triangleIndex * 9 + 4] = positions[faces[i]->positionIndices[2]-1]->y;
				vertexCoords[triangleIndex * 9 + 5] = positions[faces[i]->positionIndices[2]-1]->z;

				vertexCoords[triangleIndex * 9 + 6] = positions[faces[i]->positionIndices[3]-1]->x;
				vertexCoords[triangleIndex * 9 + 7] = positions[faces[i]->positionIndices[3]-1]->y;
				vertexCoords[triangleIndex * 9 + 8] = positions[faces[i]->positionIndices[3]-1]->z;


				vertexNormalCoords[triangleIndex * 9] =     normals[faces[i]->normalIndices[1]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 1] = normals[faces[i]->normalIndices[1]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 2] = normals[faces[i]->normalIndices[1]-1]->z;

				vertexNormalCoords[triangleIndex * 9 + 3] = normals[faces[i]->normalIndices[2]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 4] = normals[faces[i]->normalIndices[2]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 5] = normals[faces[i]->normalIndices[2]-1]->z;

				vertexNormalCoords[triangleIndex * 9 + 6] = normals[faces[i]->normalIndices[3]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 7] = normals[faces[i]->normalIndices[3]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 8] = normals[faces[i]->normalIndices[3]-1]->z;

				triangleIndex++;
			}
			else 
			{
				vertexTexCoords[triangleIndex * 6] =     texcoords[faces[i]->texcoordIndices[0]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 1] = 1-texcoords[faces[i]->texcoordIndices[0]-1]->y;
				
				vertexTexCoords[triangleIndex * 6 + 2] = texcoords[faces[i]->texcoordIndices[1]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 3] = 1-texcoords[faces[i]->texcoordIndices[1]-1]->y;

				vertexTexCoords[triangleIndex * 6 + 4] = texcoords[faces[i]->texcoordIndices[2]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 5] = 1-texcoords[faces[i]->texcoordIndices[2]-1]->y;

				vertexCoords[triangleIndex * 9] =     positions[faces[i]->positionIndices[0]-1]->x;
				vertexCoords[triangleIndex * 9 + 1] = positions[faces[i]->positionIndices[0]-1]->y;
				vertexCoords[triangleIndex * 9 + 2] = positions[faces[i]->positionIndices[0]-1]->z;

				vertexCoords[triangleIndex * 9 + 3] = positions[faces[i]->positionIndices[1]-1]->x;
				vertexCoords[triangleIndex * 9 + 4] = positions[faces[i]->positionIndices[1]-1]->y;
				vertexCoords[triangleIndex * 9 + 5] = positions[faces[i]->positionIndices[1]-1]->z;

				vertexCoords[triangleIndex * 9 + 6] = positions[faces[i]->positionIndices[2]-1]->x;
				vertexCoords[triangleIndex * 9 + 7] = positions[faces[i]->positionIndices[2]-1]->y;
				vertexCoords[triangleIndex * 9 + 8] = positions[faces[i]->positionIndices[2]-1]->z;


				vertexNormalCoords[triangleIndex * 9] =     normals[faces[i]->normalIndices[0]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 1] = normals[faces[i]->normalIndices[0]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 2] = normals[faces[i]->normalIndices[0]-1]->z;

				vertexNormalCoords[triangleIndex * 9 + 3] = normals[faces[i]->normalIndices[1]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 4] = normals[faces[i]->normalIndices[1]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 5] = normals[faces[i]->normalIndices[1]-1]->z;

				vertexNormalCoords[triangleIndex * 9 + 6] = normals[faces[i]->normalIndices[2]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 7] = normals[faces[i]->normalIndices[2]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 8] = normals[faces[i]->normalIndices[2]-1]->z;

				triangleIndex++;
			}
		}
	}

	glBindVertexArray(vao);		

	unsigned int vbo[3];		
	glGenBuffers(3, &vbo[0]);	

	glBindBuffer(GL_ARRAY_BUFFER, vbo[0]); 
	glBufferData(GL_ARRAY_BUFFER, nTriangles * 9 * sizeof(float), vertexCoords, GL_STATIC_DRAW);	    
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);     

	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]); 
	glBufferData(GL_ARRAY_BUFFER, nTriangles * 6 * sizeof(float), vertexTexCoords, GL_STATIC_DRAW);	
	glEnableVertexAttribArray(1);  
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, NULL); 
	
	glBindBuffer(GL_ARRAY_BUFFER, vbo[2]); 
	glBufferData(GL_ARRAY_BUFFER, nTriangles * 9 * sizeof(float), vertexNormalCoords, GL_STATIC_DRAW);	    
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, NULL);     
	
	delete vertexCoords;
	delete vertexTexCoords;
	delete vertexNormalCoords;
}


void PolygonalMesh::Draw()
{
	glEnable(GL_DEPTH_TEST);
	glBindVertexArray(vao); 
	glDrawArrays(GL_TRIANGLES, 0, nTriangles * 3);	
	glDisable(GL_DEPTH_TEST);
}


PolygonalMesh::~PolygonalMesh()
{
	for(unsigned int i = 0; i < rows.size(); i++) delete rows[i];   
	for(unsigned int i = 0; i < positions.size(); i++) delete positions[i];
	for(unsigned int i = 0; i < submeshFaces.size(); i++)
		for(unsigned int j = 0; j < submeshFaces.at(i).size(); j++)
			delete submeshFaces.at(i).at(j); 
	for(unsigned int i = 0; i < normals.size(); i++) delete normals[i];
	for(unsigned int i = 0; i < texcoords.size(); i++) delete texcoords[i];
}



class Shader
{
protected:
	unsigned int shaderProgram;

public:
	Shader()
	{
		shaderProgram = 0;
	}

	~Shader()
	{
		if(shaderProgram) glDeleteProgram(shaderProgram);
	}

	void Run()
	{
		if(shaderProgram) glUseProgram(shaderProgram);
	}

	virtual void UploadInvM(mat4& InVM) { }

	virtual void UploadMVP(mat4& MVP) { }
    
    virtual void UploadVP(mat4& VP) {}

	virtual void UploadSamplerID() { }
    
    virtual void UploadAttributes(vec3 ka, vec3 kd, vec3 ks, float shininess) { }
    
    virtual void UploadLightAttributes(vec3 La, vec3 Le, vec4 worldLightPosition) { }
    
    virtual void UploadSpotLightAttributes(vec3 La, vec3 Le, vec4 SpotLightPosition) { }
    
    virtual void UploadM(mat4& M) { }
    
    virtual void UploadEyePosition(vec3 wEye) { }
};



class MeshShader : public Shader
{
public:
	MeshShader()
	{
		const char *vertexSource = "\n\
			#version 410 \n\
    		precision highp float; \n\
			\n\
			in vec3 vertexPosition; \n\
			in vec2 vertexTexCoord; \n\
			in vec3 vertexNormal; \n\
			uniform mat4 M, InvM, MVP; \n\
            uniform vec3 worldEyePosition; \n\
            uniform vec4 worldLightPosition;\n\
            uniform vec4 SpotLightPosition;\n\
			out vec2 texCoord; \n\
            out vec3 worldNormal;\n\
            out vec3 worldView; \n\
            out vec3 worldLight; \n\
            out vec3 spotLight; \n\
            out vec4 pos; \n\
            out vec3 spotpos;\n\
			\n\
			void main() { \n\
				texCoord = vertexTexCoord; \n\
                vec4 worldPosition = vec4(vertexPosition, 1) * M;\n\
                worldLight = worldLightPosition.xyz * worldPosition.w - worldPosition.xyz * worldLightPosition.w; \n\
                spotLight = SpotLightPosition.xyz * worldPosition.w - worldPosition.xyz * SpotLightPosition.w; \n\
                worldView = worldEyePosition - worldPosition.xyz; \n\
				worldNormal = (InvM * vec4(vertexNormal, 0.0)).xyz; \n\
                pos = worldPosition; \n\
                spotpos = spotLight; \n\
				gl_Position = vec4(vertexPosition, 1) * MVP; \n\
			} \n\
		"; 

		const char *fragmentSource = "\n\
			#version 410 \n\
    		precision highp float; \n\
			\n\
			uniform sampler2D samplerUnit; \n\
            uniform vec3 La, Le;\n\
            uniform vec3 SpotLa, SpotLe;\n\
            uniform vec3 ka, kd, ks;\n\
            uniform float shininess;\n\
			in vec2 texCoord; \n\
            in vec3 worldNormal;\n\
            in vec3 worldView;\n\
            in vec3 worldLight;\n\
            in vec3 spotLight; \n\
            in vec4 pos; \n\
            in vec3 spotpos; \n\
			out vec4 fragmentColor; \n\
			\n\
			void main() { \n\
                vec3 N = normalize(worldNormal); \n\
                vec3 V = normalize(worldView); \n\
				vec3 L = normalize(worldLight); \n\
                vec3 H = normalize(V + L); \n\
                vec3 Ls = normalize(spotLight); \n\
                vec3 Hs = normalize(V + Ls); \n\
                vec3 texel = texture(samplerUnit, texCoord).xyz; \n\
                vec3 color = La * ka + Le * kd * texel * max(0.0, dot(L,N)) + Le * ks * pow(max(0.0, dot(H,N)), shininess);\n\
                vec3 spotcolor = 4*(SpotLa * ka + SpotLe * kd * texel * max(0.0, dot(Ls,N)) + SpotLe * ks * pow(max(0.0, dot(Hs,N)), shininess))/dot(pos.xyz-spotpos,pos.xyz-spotpos);\n\
                fragmentColor = vec4(color+spotcolor,1); \n\
			} \n\
		";

		unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
		if (!vertexShader) { printf("Error in vertex shader creation\n"); exit(1); }

		glShaderSource(vertexShader, 1, &vertexSource, NULL);
		glCompileShader(vertexShader);
		checkShader(vertexShader, "Vertex shader error");

		unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
		if (!fragmentShader) { printf("Error in fragment shader creation\n"); exit(1); }

		glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
		glCompileShader(fragmentShader);
		checkShader(fragmentShader, "Fragment shader error");

		shaderProgram = glCreateProgram();
		if (!shaderProgram) { printf("Error in shader program creation\n"); exit(1); }

		glAttachShader(shaderProgram, vertexShader);
		glAttachShader(shaderProgram, fragmentShader);

		glBindAttribLocation(shaderProgram, 0, "vertexPosition");
		glBindAttribLocation(shaderProgram, 1, "vertexTexCoord");
		glBindAttribLocation(shaderProgram, 2, "vertexNormal");

		glBindFragDataLocation(shaderProgram, 0, "fragmentColor");

		glLinkProgram(shaderProgram);
		checkLinking(shaderProgram);
	}

	void UploadSamplerID()
	{
		int samplerUnit = 0; 
		int location = glGetUniformLocation(shaderProgram, "samplerUnit");
		glUniform1i(location, samplerUnit);
		glActiveTexture(GL_TEXTURE0 + samplerUnit); 
	}

	void UploadInvM(mat4& InvM)
	{
		int location = glGetUniformLocation(shaderProgram, "InvM");
		if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, InvM); 
		else printf("uniform InvM cannot be set\n");
	}

	void UploadMVP(mat4& MVP)
	{
		int location = glGetUniformLocation(shaderProgram, "MVP");
		if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, MVP);
		else printf("uniform MVP cannot be set\n");
	}
    
    void UploadM(mat4& M)
    {
        int location = glGetUniformLocation(shaderProgram, "M");
        if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, M);
        else printf("uniform M cannot be set\n");
    }
    
    // Upload attributes from Material class
    void UploadAttributes(vec3 ka, vec3 kd, vec3 ks, float shininess){
        int location = glGetUniformLocation(shaderProgram, "ka");
        if (location >= 0) glUniform3f(location, ka.x, ka.y, ka.z);
        else printf("unifrom ka cannot be set\n");
        
        location = glGetUniformLocation(shaderProgram, "kd");
        if (location >= 0) glUniform3f(location, kd.x, kd.y, kd.z);
        else printf("unifrom kd cannot be set\n");
        
        location = glGetUniformLocation(shaderProgram, "ks");
        if (location >= 0) glUniform3f(location, ks.x, ks.y, ks.z);
        else printf("unifrom ks cannot be set\n");
        
        location = glGetUniformLocation(shaderProgram, "shininess");
        if (location >= 0) glUniform1f(location, shininess);
        else printf("uniform shininess cannot be set\n");
    }
    
    // Upload attributes from Light class
    void UploadLightAttributes(vec3 La, vec3 Le, vec4 worldLightPosition){
        int location = glGetUniformLocation(shaderProgram, "La");
        if( location >=0) glUniform3f(location, La.x, La.y, La.z);
        else printf("unifrom La cannot be set\n");
        
        location = glGetUniformLocation(shaderProgram, "Le");
        if( location >=0) glUniform3f(location, Le.x, Le.y, Le.z);
        else printf("unifrom Le cannot be set\n");
        
        location = glGetUniformLocation(shaderProgram, "worldLightPosition");
        if( location >=0) glUniform4fv(location, 1, &worldLightPosition.v[0]);
        else printf("unifrom worldLightPosition cannot be set\n");
    }
    
    // Upload attributes from SpotLight class
    void UploadSpotLightAttributes(vec3 La, vec3 Le, vec4 SpotLightPosition){
        int location = glGetUniformLocation(shaderProgram, "SpotLa");
        if( location >=0) glUniform3f(location, La.x, La.y, La.z);
        else printf("unifrom spotLa cannot be set\n");
        
        location = glGetUniformLocation(shaderProgram, "SpotLe");
        if( location >=0) glUniform3f(location, Le.x, Le.y, Le.z);
        else printf("unifrom spotLe cannot be set\n");
        
        location = glGetUniformLocation(shaderProgram, "SpotLightPosition");
        if( location >=0) glUniform4fv(location, 1, &SpotLightPosition.v[0]);
        else printf("unifrom SpotLightPosition cannot be set\n");
    }
    
    // Upload attributes from Camera class
    void UploadEyePosition(vec3 wEye){
        int location = glGetUniformLocation(shaderProgram, "worldEyePosition");
        if(location >=0) glUniform3f(location, wEye.x, wEye.y, wEye.z);
        else printf("unifrom worldEyePosition cannot be set\n");
    }
};



class InfiniteQuadShader : public Shader
{
public:
    InfiniteQuadShader()
    {
        const char *vertexSource = "\n\
        #version 410 \n\
        precision highp float; \n\
        \n\
        in vec4 vertexPosition; \n\
        in vec2 vertexTexCoord; \n\
        in vec3 vertexNormal; \n\
        uniform mat4 M, InvM, MVP; \n\
        out vec2 texCoord; \n\
        out vec4 worldPosition;\n\
        out vec3 worldNormal;\n\
        \n\
        void main() { \n\
            texCoord = vertexTexCoord; \n\
            worldPosition = vertexPosition * M;\n\
            worldNormal = (InvM * vec4(vertexNormal, 0.0)).xyz; \n\
            gl_Position = vertexPosition * MVP; \n\
        } \n\
        ";
        
        const char *fragmentSource = "\n\
        #version 410 \n\
        precision highp float; \n\
        \n\
        uniform sampler2D samplerUnit; \n\
        uniform vec3 La, Le;\n\
        uniform vec3 ka, kd, ks;\n\
        uniform float shininess;\n\
        uniform vec3 worldEyePosition;\n\
        uniform vec4 worldLightPosition;\n\
        in vec2 texCoord;\n\
        in vec4 worldPosition;\n\
        in vec3 worldNormal;\n\
        out vec4 fragmentColor; \n\
        \n\
        void main() { \n\
            vec3 N = normalize(worldNormal); \n\
            vec3 V = normalize(worldEyePosition * worldPosition.w - worldPosition.xyz); \n\
            vec3 L = normalize(worldLightPosition.xyz * worldPosition.w - worldPosition.xyz * worldLightPosition.w); \n\
            vec3 H = normalize(V + L); \n\
            vec2 position = worldPosition.xz / worldPosition.w;\n\
            vec2 tex = position.xy - floor(position.xy); \n\
            vec3 texel = texture(samplerUnit, tex).xyz; \n\
            vec3 color = La * ka + Le * kd * texel * max(0.0, dot(L,N)) + Le * ks * pow(max(0.0, dot(H,N)), shininess);\n\
            fragmentColor = vec4(color,1); \n\
        } \n\
        ";
        
        unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
        if (!vertexShader) { printf("Error in vertex shader creation\n"); exit(1); }
        
        glShaderSource(vertexShader, 1, &vertexSource, NULL);
        glCompileShader(vertexShader);
        checkShader(vertexShader, "Vertex shader error");
        
        unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        if (!fragmentShader) { printf("Error in fragment shader creation\n"); exit(1); }
        
        glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
        glCompileShader(fragmentShader);
        checkShader(fragmentShader, "Fragment shader error");
        
        shaderProgram = glCreateProgram();
        if (!shaderProgram) { printf("Error in shader program creation\n"); exit(1); }
        
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        
        glBindAttribLocation(shaderProgram, 0, "vertexPosition");
        glBindAttribLocation(shaderProgram, 1, "vertexTexCoord");
        glBindAttribLocation(shaderProgram, 2, "vertexNormal");
        
        glBindFragDataLocation(shaderProgram, 0, "fragmentColor");
        
        glLinkProgram(shaderProgram);
        checkLinking(shaderProgram);
    }
    
    void UploadSamplerID()
    {
        int samplerUnit = 0;
        int location = glGetUniformLocation(shaderProgram, "samplerUnit");
        glUniform1i(location, samplerUnit);
        glActiveTexture(GL_TEXTURE0 + samplerUnit);
    }
    
    void UploadInvM(mat4& InvM)
    {
        int location = glGetUniformLocation(shaderProgram, "InvM");
        if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, InvM);
        else printf("uniform InvM cannot be set\n");
    }
    
    void UploadMVP(mat4& MVP)
    {
        int location = glGetUniformLocation(shaderProgram, "MVP");
        if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, MVP);
        else printf("uniform MVP cannot be set\n");
    }
    
    void UploadM(mat4& M)
    {
        int location = glGetUniformLocation(shaderProgram, "M");
        if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, M);
        else printf("uniform M cannot be set\n");
    }
    
    // Upload attributes from Material class
    void UploadAttributes(vec3 ka, vec3 kd, vec3 ks, float shininess){
        int location = glGetUniformLocation(shaderProgram, "ka");
        if (location >= 0) glUniform3f(location, ka.x, ka.y, ka.z);
        else printf("unifrom ka cannot be set\n");
        
        location = glGetUniformLocation(shaderProgram, "kd");
        if (location >= 0) glUniform3f(location, kd.x, kd.y, kd.z);
        else printf("unifrom kd cannot be set\n");
        
        location = glGetUniformLocation(shaderProgram, "ks");
        if (location >= 0) glUniform3f(location, ks.x, ks.y, ks.z);
        else printf("unifrom ks cannot be set\n");
        
        location = glGetUniformLocation(shaderProgram, "shininess");
        if (location >= 0) glUniform1f(location, shininess);
        else printf("uniform shininess cannot be set\n");
    }
    
    // Upload attributes from Light class
    void UploadLightAttributes(vec3 La, vec3 Le, vec4 worldLightPosition){
        int location = glGetUniformLocation(shaderProgram, "La");
        if( location >=0) glUniform3f(location, La.x, La.y, La.z);
        else printf("unifrom La cannot be set\n");
        
        location = glGetUniformLocation(shaderProgram, "Le");
        if( location >=0) glUniform3f(location, Le.x, Le.y, Le.z);
        else printf("unifrom Le cannot be set\n");
        
        location = glGetUniformLocation(shaderProgram, "worldLightPosition");
        if( location >=0) glUniform4fv(location, 1, &worldLightPosition.v[0]);
        else printf("unifrom worldLightPosition cannot be set\n");
    }
    
    // Upload attributes from Camera class
    void UploadEyePosition(vec3 wEye){
        int location = glGetUniformLocation(shaderProgram, "worldEyePosition");
        if(location >=0) glUniform3f(location, wEye.x, wEye.y, wEye.z);
        else printf("unifrom worldEyePosition cannot be set\n");
    }
};


class ShadowShader : public Shader
{
public:
    ShadowShader()
    {
        const char *vertexSource = "\n\
        #version 410 \n\
        precision highp float; \n\
        \n\
        in vec3 vertexPosition; \n\
        in vec2 vertexTexCoord; \n\
        in vec3 vertexNormal; \n\
        uniform mat4 M, VP; \n\
        uniform vec4 worldLightPosition;\n\
        void main() { \n\
            vec4 p = vec4(vertexPosition,1) * M; \n\
            vec3 s; \n\
            s.y = -0.999; \n\
            s.x = (p.x - worldLightPosition.x) / (p.y - worldLightPosition.y) * (s.y - worldLightPosition.y) + worldLightPosition.x; \n\
            s.z = (p.z - worldLightPosition.z) / (p.y - worldLightPosition.y) * (s.y - worldLightPosition.y) + worldLightPosition.z; \n\
            gl_Position = vec4(s, 1) * VP; \n\
        } \n\
        ";
        
        const char *fragmentSource = " \n\
        #version 410 \n\
        precision highp float; \n\
        \n\
        out vec4 fragmentColor; \n\
        \n\
        void main() { \n\
            fragmentColor = vec4(0.0, 0.1, 0.0, 1); \n\
        } \n\
        ";
        
        unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
        if (!vertexShader) { printf("Error in vertex shader creation\n"); exit(1); }
        
        glShaderSource(vertexShader, 1, &vertexSource, NULL);
        glCompileShader(vertexShader);
        checkShader(vertexShader, "Vertex shader error");
        
        unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        if (!fragmentShader) { printf("Error in fragment shader creation\n"); exit(1); }
        
        glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
        glCompileShader(fragmentShader);
        checkShader(fragmentShader, "Fragment shader error");
        
        shaderProgram = glCreateProgram();
        if (!shaderProgram) { printf("Error in shader program creation\n"); exit(1); }
        
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        
        glBindAttribLocation(shaderProgram, 0, "vertexPosition");
        glBindAttribLocation(shaderProgram, 1, "vertexTexCoord");
        glBindAttribLocation(shaderProgram, 2, "vertexNormal");
        
        glBindFragDataLocation(shaderProgram, 0, "fragmentColor");
        
        glLinkProgram(shaderProgram);
        checkLinking(shaderProgram);
    }
    
    
    void UploadVP(mat4& VP)
    {
        int location = glGetUniformLocation(shaderProgram, "VP");
        if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, VP);
        else printf("uniform VP cannot be set\n");
    }
    
    void UploadM(mat4& M)
    {
        int location = glGetUniformLocation(shaderProgram, "M");
        if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, M);
        else printf("uniform M cannot be set\n");
    }
    
    
    // Upload attributes from Light class
    void UploadLightAttributes(vec3 La, vec3 Le, vec4 worldLightPosition){
        
        int location = glGetUniformLocation(shaderProgram, "worldLightPosition");
        if( location >=0) glUniform4fv(location, 1, &worldLightPosition.v[0]);
        else printf("unifrom worldLightPosition cannot be set\n");
    }
    
};



extern "C" unsigned char* stbi_load(char const *filename, int *x, int *y, int *comp, int req_comp);

class Texture
{
	unsigned int textureId;

public:
	Texture(const std::string& inputFileName)
	{
		unsigned char* data;
		int width; int height; int nComponents = 4;
		
		data = stbi_load(inputFileName.c_str(), &width, &height, &nComponents, 0);

		if(data == NULL) 
		{ 
			return;
		}

		glGenTextures(1, &textureId); 
		glBindTexture(GL_TEXTURE_2D, textureId); 
		
		if(nComponents == 3) glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
		if(nComponents == 4) glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

		delete data; 
	}

	void Bind()
	{
		glBindTexture(GL_TEXTURE_2D, textureId);
	}
};



class Light
{
    vec3 La, Le;
    vec4 worldLightPosition;
    
public:
    Light(vec3 inLa, vec3 inLe, vec4 inWorldLightPosition){
        La = inLa;
        Le = inLe;
        worldLightPosition = inWorldLightPosition;
    }
    
    void UploadAttributes(Shader* s){
        s->UploadLightAttributes(La, Le, worldLightPosition);
    }
    
    void SetPointLightSource(vec3& pos){
        worldLightPosition.v[0] = pos.x;
        worldLightPosition.v[1] = pos.y;
        worldLightPosition.v[2] = pos.z;
        worldLightPosition.v[3] = 1.0;
    }
    
    void SetDirectionalLightSource(vec3& dir){
        worldLightPosition.v[0] = dir.x;
        worldLightPosition.v[1] = dir.y;
        worldLightPosition.v[2] = dir.z;
        worldLightPosition.v[3] = 0.0;
    }
    
    
};

class SpotLight
{
    vec3 La, Le;
    vec4 worldLightPosition;
    
public:
    SpotLight(vec3 inLa, vec3 inLe, vec4 inWorldLightPosition){
        La = inLa;
        Le = inLe;
        worldLightPosition = inWorldLightPosition;
    }
    
    void UploadAttributes(Shader* s){
        s->UploadSpotLightAttributes(La, Le, worldLightPosition);
    }
    
    void Move(vec3 velocity, float dt){
        worldLightPosition = worldLightPosition + vec4(velocity.x*dt,velocity.y*dt,velocity.z*dt,0);
    }
    
    void SetPointLightSource(vec3& pos){
        worldLightPosition.v[0] = pos.x;
        worldLightPosition.v[1] = pos.y;
        worldLightPosition.v[2] = pos.z;
        worldLightPosition.v[3] = 1.0;
    }
    
    
};



Light* light;
SpotLight* spotLight;


class Material
{
	Shader* shader;
	Texture* texture;

    vec3 ka;
    vec3 kd;
    vec3 ks;
    float shininess;
    
public:
	Material(Shader* s, vec3 inka, vec3 inkd, vec3 inks, float inshininess, Texture* t = 0)
	{
		shader = s;
		texture = t;
        ka = inka;
        kd = inkd;
        ks = inks;
        shininess = inshininess;
	}

	Shader* GetShader() { return shader; }

	void UploadAttributes()
	{
		if(texture)
		{
            shader->UploadAttributes(ka, kd, ks, shininess);
			shader->UploadSamplerID();
			texture->Bind();
		}

	}
    
};

class Mesh
{
	Geometry* geometry;
	Material* material;

public:
	Mesh(Geometry* g, Material* m)
	{
		geometry = g;
		material = m;
	}

	Shader* GetShader() { return material->GetShader(); }

	void Draw()
	{
		material->UploadAttributes();
		geometry->Draw();
	}
};

vec3 initialPos; // used for heart TrackingShot
vec3 initialIcePos;

class Camera {
   vec3  wEye, wLookat, wVup;
   float fov, asp, fp, bp;
    
    vec3 v;   // moving speed of 'w' and 's'
    double angularV;    //rotating angular velocity of 'a' and 'd'

public:
	Camera()
	{
		wEye = vec3(0.0, 0.0, 2.0);
		wLookat = vec3(0.0, 0.0, 0.0);
		wVup = vec3(0.0, 1.0, 0.0);
		fov = M_PI / 4.0; asp = 1.0; fp = 0.01; bp = 10.0;		
	}
	
	void SetAspectRatio(float a) { asp = a; }

	mat4 GetViewMatrix() 
	{ 
		vec3 w = (wEye - wLookat).normalize();
		vec3 u = cross(wVup, w).normalize();
		vec3 v = cross(w, u);
	
		return  
			mat4(	
				1.0f,    0.0f,    0.0f,    0.0f,
				0.0f,    1.0f,    0.0f,    0.0f,
				0.0f,    0.0f,    1.0f,    0.0f,
				-wEye.x, -wEye.y, -wEye.z, 1.0f ) *
			mat4(	
				u.x,  v.x,  w.x,  0.0f,
				u.y,  v.y,  w.y,  0.0f,
				u.z,  v.z,  w.z,  0.0f,
				0.0f, 0.0f, 0.0f, 1.0f );
   }

	mat4 GetProjectionMatrix() 
	{ 
		float sy = 1/tan(fov/2);
		return mat4(
			sy/asp, 0.0f,  0.0f,               0.0f,
			0.0f,   sy,    0.0f,               0.0f,
			0.0f,   0.0f, -(fp+bp)/(bp - fp), -1.0f,
			0.0f,   0.0f, -2*fp*bp/(bp - fp),  0.0f);
	}
    
    void Control(){
        //the 2 speeds include directions
        float speed = 1;
        vec3 dif = wLookat - wEye;
        v = dif * (keyboardState['w'] * speed + keyboardState['s'] * -speed);
        
        float angularSpeed = 1;
        angularV = keyboardState['a'] * -angularSpeed + keyboardState['d'] * angularSpeed;
    }
    
    void Move(float dt) {
        //move the obj using v and angularV
        wEye = wEye + v * dt;
        wLookat = wLookat + v * dt;
        
        vec3 w = (wLookat - wEye).normalize();
        vec3 v = wVup.normalize();
        vec3 u = cross(w, v);
        
        vec3 wPrime = (w*cos(angularV*dt) + u * sin(angularV*dt)) * (wLookat - wEye).length();
        //wLookat = wPrime + wEye;
        wEye = wLookat - wPrime;
    }
    
    void SetwEye(vec3 pos){
        wEye = pos;
    }
    
    vec3 GetwEye(){
        return wEye;
    }
    
    vec3 GetLookat(){
        return wLookat;
    }
    
    void UploadAttributes(Shader* s){
        s->UploadEyePosition(wEye);
    }
    
    void TrackingShot(float t)
    {
        wEye = heart(t);
    }
    
    vec3 heart(float t){
        vec3 v0 = vec3(0,0,-2).normalize();
        
        vec3 v1;
        if(keyboardState['t']){
            v1 = (wLookat-initialPos).normalize();
        }else{
            v1 = (wLookat-wEye).normalize();
        }
        
        float co = v0.z*v1.z; // dot product
        float si = cross(v0, v1).length(); // size of cross product
        
        float heartToEyeX = x(t);
        float heartToEyeZ = z(t)-2;
        
        vec3 eyeToPos = vec3(heartToEyeX*co-heartToEyeZ*si, 0 ,heartToEyeX*si+heartToEyeZ*co);
        
        return eyeToPos+initialPos;
    }
    
    
    //parametric function scaled down by 3
    float x(float t){
        return 16*sin(t)*sin(t)*sin(t)/3;
    }
    
    float z(float t){
        t += M_PI;
        return -( 13*cos(t)-5*cos(2*t)-2*cos(3*t)-cos(4*t) )/3;
    }
};

Camera camera;

class Object
{
	Shader* shader;
	Mesh *mesh;

	vec3 position;
	vec3 scaling;
    
	float orientation;
    float angularV;
    vec3 velocity = vec3(0,0,0);
    
    float rotation;

public:
	Object(Mesh *m, vec3 position = vec3(0.0, 0.0, 0.0), vec3 scaling = vec3(1.0, 1.0, 1.0), float orientation = 0.0, float rotation = 0.0) : position(position), scaling(scaling), orientation(orientation), rotation(rotation)
	{
		shader = m->GetShader();
		mesh = m;
	}

	vec3 GetPosition() { return position; }
    
    vec3 GetVelocity(){return velocity;}
    
    void MovePosition(float dt){
        
        vec3 dif = camera.GetLookat()-camera.GetwEye();
        velocity = dif * (keyboardState['w'] - keyboardState['s']);
        position = position + velocity * dt;
        
        if(keyboardState['a'] or keyboardState['d']){
            float angularSpeed = 180/M_PI;
            angularV = keyboardState['a'] * -angularSpeed + keyboardState['d'] * angularSpeed;
            orientation += angularV * dt;
            
        }

    }
    
    void Snake(float t, float dt){ //x=R*cos(t), y=A*sin(3*t), z=R*sin(t), R=3, A=0.5
        velocity = vec3(-3*sin(t), 1.5*cos(3*t), 3*cos(t));
        position = position + velocity * dt;
    }
    
    void frenet(float dt){
        if(keyboardState['d']){
            if(keyboardState['w']){
                rotation -= 60*dt;
                rotation = fmax(rotation, -45);

            } else if (keyboardState['s']){
                rotation += 60*dt;
                rotation = fmin(rotation, 45);

            } else {
                if(rotation>0){
                    rotation -= 60*dt;
                    rotation = fmax(rotation, 0);

                } else {
                    rotation += 60*dt;
                    rotation = fmin(rotation, 0);

                }
            }
        } else if(keyboardState['a']){
            if(keyboardState['w']){
                rotation += 60*dt;
                rotation = fmin(rotation, 45);
                
            }else if(keyboardState['s']){
                rotation -= 60*dt;
                rotation = fmax(rotation, -45);

            } else {
                if(rotation>0){
                    rotation -= 60*dt;
                    rotation = fmax(rotation, 0);

                } else {
                    rotation += 60*dt;
                    rotation = fmin(rotation, 0);

                }
            }
            
        } else {
            if(rotation>0){
                rotation -= 60*dt;
                rotation = fmax(rotation, 0);
            } else {
                rotation += 60*dt;
                rotation = fmin(rotation, 0);

            }
        }
        
    }
    
    void Spin(float dt){
        orientation += dt*200;
    }
    
    void Rotate(float dt){
        rotation += dt*200;
    }

	void Draw()
	{
		shader->Run();
        spotLight->UploadAttributes(shader);
        light->UploadAttributes(shader);
        
		UploadAttributes(shader);
		mesh->Draw();
	}
    
    void DrawShadow(Shader* shadowShader){
        shadowShader->Run();
        
        UploadAttributes(shadowShader);
        
        light->UploadAttributes(shadowShader);
        //spotLight->UploadAttributes(shadowShader);
        camera.UploadAttributes(shadowShader);
        
        mesh->Draw();
    }

	void UploadAttributes(Shader* shader)
	{
		mat4 T = mat4(
			1.0,			0.0,			0.0,			0.0,
			0.0,			1.0,			0.0,			0.0,
			0.0,			0.0,			1.0,			0.0,
			position.x,		position.y,		position.z,		1.0);

		mat4 InvT = mat4(
			1.0,			0.0,			0.0,			0.0,
			0.0,			1.0,			0.0,			0.0,
			0.0,			0.0,			1.0,			0.0,
			-position.x,	-position.y,	-position.z,	1.0);

		mat4 S = mat4(
			scaling.x,		0.0,			0.0,			0.0,
			0.0,			scaling.y,		0.0,			0.0,
			0.0,			0.0,			scaling.z,		0.0,
			0.0,			0.0,			0.0,			1.0);

		mat4 InvS = mat4(
			1.0/scaling.x,	0.0,			0.0,			0.0,
			0.0,			1.0/scaling.y,	0.0,			0.0,
			0.0,			0.0,			1.0/scaling.z,	0.0,
			0.0,			0.0,			0.0,			1.0);

		float alpha = orientation / 180.0 * M_PI;
        float beta = rotation / 180.0 * M_PI;
        
		mat4 R = mat4(
			cos(alpha),		0.0,			sin(alpha),		0.0,
			0.0,			1.0,			0.0,			0.0,
			-sin(alpha),	0.0,			cos(alpha),		0.0,
			0.0,			0.0,			0.0,			1.0)

        ;

		mat4 InvR =
        mat4(
			cos(alpha),		0.0,			-sin(alpha),	0.0,
			0.0,			1.0,			0.0,			0.0,
			sin(alpha),		0.0,			cos(alpha),		0.0,
			0.0,			0.0,			0.0,			1.0);
        
        
        vec3 u;
        if(not keyboardState['t']){
            u = (camera.GetLookat()-camera.GetwEye()).normalize(); //axis from which the obj rotates u.y=0
        } else {
            u = (camera.GetLookat()-initialPos).normalize(); //axis from which the obj rotates u.y=0
        }
        
        
        mat4 Rz =
        mat4(cos(beta) + u.x*u.x*(1-cos(beta)) ,  -u.z*sin(beta),    u.x*u.z*(1-cos(beta)),     0,
             u.z*sin(beta),                       cos(beta),         -u.x*sin(beta),            0,
             u.z*u.x*(1-cos(beta)),               u.x*sin(beta),     cos(beta)+u.z*u.z*(1-cos(beta)), 0,
             0, 0 , 0 , 1
        );
        
        mat4 InvRz =
        mat4(cos(beta) + u.x*u.x*(1-cos(beta)) ,  u.z*sin(beta) ,    u.z*u.x*(1-cos(beta)) ,     0,
             -u.z*sin(beta),                       cos(beta),         u.x*sin(beta),            0,
             u.x*u.z*(1-cos(beta)),               -u.x*sin(beta),     cos(beta)+u.z*u.z*(1-cos(beta)), 0,
             0, 0 , 0 , 1
             );


		mat4 M = S * R * Rz * T;
		mat4 InvM = InvT * InvRz * InvR * InvS;

		mat4 MVP = M * camera.GetViewMatrix() * camera.GetProjectionMatrix();
        mat4 VP = camera.GetViewMatrix() * camera.GetProjectionMatrix();
        
        shader->UploadVP(VP);
		shader->UploadInvM(InvM);
		shader->UploadMVP(MVP);
        shader->UploadM(M);
	}
};



Object* thunderbolt;
Object* airscrew;
Object* heli;
Object* rotor;
Object* heli1;
Object* rotor1;
Object* heli2;
Object* rotor2;

class Scene
{
	MeshShader *meshShader;
    InfiniteQuadShader *infShader;
    ShadowShader *shadowShader;
	
	std::vector<Texture*> textures;
	std::vector<Material*> materials;
	std::vector<Geometry*> geometries;
	std::vector<Mesh*> meshes;
	std::vector<Object*> objects;

public:
	Scene() 
	{ 
		meshShader = 0;
	}

	void Initialize()
	{
        //5 10 5 0
        light = new Light(vec3(1, 1, 1), vec3(1, 1, 1), vec4(5.0, 10.0, 5.0, 0.0)); //directional light source
        spotLight = new SpotLight(vec3(1.0, 1.0, 1.0), vec3(1.0, 1.0, 1.0), vec4(0.0, 2.0, 0.0, 1.0)); //point light source above the avatar
        
		meshShader = new MeshShader();
        infShader = new InfiniteQuadShader();
        shadowShader = new ShadowShader();
        
        //wEye is at 0,0,2
        //                  0             1            2           3          4        5       6
        // texture:        tb        tree       airscrew        heli       rotor    lava
        // material:       tb     smalltrees    airscrew      bigtree     heli    rotor   lava
        // geometries:    tb.obj   tree.obj     airscrew       heli          rotor    infiniteQuad
        // meshes:         tb     smalltrees    airscrew      bigtree     heli    quad
        // object:         tb     smalltrees    airscrew       heli        lava
        
        //thunderbolt obj
		textures.push_back(new Texture("orangeblue.png"));
		materials.push_back(new Material(meshShader, vec3(0.1,0.1,0.1), vec3(0.9,0.9,0.9), vec3(0.0,0.0,0.0), 0, textures[0]));
		geometries.push_back(new PolygonalMesh("thunderbolt_body.obj"));
		meshes.push_back(new Mesh(geometries[0], materials[0]));
		
		thunderbolt = new Object(meshes[0], vec3(0.0, -0.5, 0), vec3(0.02, 0.02, 0.02), 0.0, 0.0);
		objects.push_back(thunderbolt);
        
        
        //3 tree objs diffuse
        textures.push_back(new Texture("tree.png"));
        materials.push_back(new Material(meshShader, vec3(0.1,0.1,0.1), vec3(0.9,0.9,0.9), vec3(0.0,0.0,0.0), 0, textures[1]));
        geometries.push_back(new PolygonalMesh("tree.obj"));
        meshes.push_back(new Mesh(geometries[1], materials[1]));
        
        float depth = -1;

        objects.push_back(new Object(meshes[1], vec3(1, depth, -1), vec3(0.01, 0.01, 0.01), 30.0));
        objects.push_back(new Object(meshes[1], vec3(2, depth, -1.5), vec3(0.01, 0.01, 0.01), 80.0));
        objects.push_back(new Object(meshes[1], vec3(-1, depth, -0.5), vec3(0.01, 0.01, 0.01), 105.0));
        
        
        
        //rotating airscrew
        textures.push_back(new Texture("orange.png"));
        materials.push_back(new Material(meshShader, vec3(0.1,0.1,0.1), vec3(0.6,0.6,0.6), vec3(0.3,0.3,0.3), 50, textures[2]));
        geometries.push_back(new PolygonalMesh("thunderbolt_airscrew.obj"));
        meshes.push_back(new Mesh(geometries[2], materials[2]));
        airscrew = new Object(meshes[2], vec3(0.0, -0.5, 0.0), vec3(0.02, 0.02, 0.02), 0.0);
        objects.push_back(airscrew);
        
        
        //big specular tree
        materials.push_back(new Material(meshShader, vec3(0.1,0.1,0.1), vec3(0.6,0.6,0.6), vec3(0.3,0.3,0.3), 50, textures[1]));
        meshes.push_back(new Mesh(geometries[1], materials[3]));
        objects.push_back(new Object(meshes[3], vec3(-0.5, depth, -3), vec3(0.1, 0.1, 0.1), 15));
        
        
        //helicoptor and rotor
        textures.push_back(new Texture("heliait.png"));
        materials.push_back(new Material(meshShader, vec3(0.1,0.1,0.1), vec3(0.6,0.6,0.6), vec3(0.3,0.3,0.3), 50, textures[3]));
        geometries.push_back(new PolygonalMesh("heli.obj"));
        meshes.push_back(new Mesh(geometries[3], materials[4]));
        
        heli = new Object(meshes[4], vec3(1, 0.5, -1), vec3(0.02, 0.02, 0.02), 45.0);
        heli1 = new Object(meshes[4], vec3(1, 0.5, -1), vec3(0.02, 0.02, 0.02), 60.0);
        heli2 = new Object(meshes[4], vec3(1, 0.5, -1), vec3(0.02, 0.02, 0.02), 30.0);
        objects.push_back(heli);
        objects.push_back(heli1);
        objects.push_back(heli2);
        
        textures.push_back(new Texture("rotor.png"));
        materials.push_back(new Material(meshShader, vec3(0.1,0.1,0.1), vec3(0.6,0.6,0.6), vec3(0.3,0.3,0.3), 50, textures[4]));
        geometries.push_back(new PolygonalMesh("mainrotor.obj"));
        meshes.push_back(new Mesh(geometries[4], materials[3]));
        
        rotor = new Object(meshes[5], vec3(0.9, 0.8, -1), vec3(0.02, 0.02, 0.02), 90.0);
//        rotor1 = new Object(meshes[5], vec3(1.9, 0.8, 1), vec3(0.02, 0.02, 0.02), 90.0);
//        rotor2 = new Object(meshes[5], vec3(2.9, 0.8, 1.5), vec3(0.02, 0.02, 0.02), 90.0);
//        objects.push_back(rotor);
//        objects.push_back(rotor1);
//        objects.push_back(rotor2);
        
        
        
        //ground as InfiniteTexturedQuad
        textures.push_back(new Texture("lava.png"));
        materials.push_back(new Material(infShader, vec3(0.1,0.1,0.1), vec3(0.9,0.9,0.9), vec3(0.0,0.0,0.0), 0, textures[5]));
        geometries.push_back(new InfiniteTexturedQuad());
        meshes.push_back(new Mesh(geometries[5], materials[6]));
        objects.push_back(new Object(meshes[6], vec3(0, -1, 0), vec3(1, 1, 1), 0));
        
        
        
	}


	~Scene()
	{
		for(int i = 0; i < textures.size(); i++) delete textures[i];
		for(int i = 0; i < materials.size(); i++) delete materials[i];
		for(int i = 0; i < geometries.size(); i++) delete geometries[i];
		for(int i = 0; i < meshes.size(); i++) delete meshes[i];
		for(int i = 0; i < objects.size(); i++) delete objects[i];
		
		if(meshShader) delete meshShader;
	}

	void Draw()
	{
		for(int i = 0; i < objects.size(); i++)
        {
            objects[i]->Draw();
            if(i != objects.size() - 1 ){
                objects[i]->DrawShadow(shadowShader);
            }
        }
	}
};

Scene scene;

void onInitialization() 
{
	glViewport(0, 0, windowWidth, windowHeight);

	scene.Initialize();
}

void onExit() 
{
	printf("exit");
}

void onDisplay() 
{	
	
	glClearColor(0, 0, 1.0, 0); 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
	
	scene.Draw();
	glutSwapBuffers(); 
	
}

bool tPressed = false;
double tHeart = 0;

void onKeyboard(unsigned char key, int x, int y)
{
	keyboardState[key] = true;
    
    if(key =='a'){
        keyboardState['d']=false;
    } else if (key == 'd'){
        keyboardState['a']=false;
    }
    
    if(key == 't' and not tPressed){
        tPressed = true;
        initialPos = camera.GetwEye();
    }
}

void onKeyboardUp(unsigned char key, int x, int y)
{

    keyboardState[key] = false;
    
    if(key == 't'){
        tPressed = false;
        camera.SetwEye(initialPos);
        tHeart = 0;
    }

}

void onReshape(int winWidth, int winHeight) 
{
	camera.SetAspectRatio((float)winWidth / winHeight);
	glViewport(0, 0, winWidth, winHeight);
}


void onIdle( ) {
    
    double t = glutGet(GLUT_ELAPSED_TIME) * 0.001;
    static double lastTime = 0.0;
    double dt = t - lastTime;
    lastTime = t;
    
    heli->Snake(t,dt);
    heli1->Snake(t+1,dt);
    heli2->Snake(t+2,dt);
    

    if(keyboardState['t'] and tPressed){
        tHeart+=dt;
        camera.TrackingShot(tHeart);
    }
    
//    vec3 pos = thunderbolt->GetPosition();
//    spotLight = new SpotLight(vec3(1.0, 1.0, 1.0), vec3(1.0, 1.0, 1.0), vec4(pos.x, 0.5, pos.z, 1.0));
    spotLight->Move(thunderbolt->GetVelocity(), dt);
    
//    rotor->Spin(2*dt);
//    rotor1->Spin(2*dt);
//    rotor2->Spin(2*dt);
    airscrew->Rotate(3*dt);
    
    camera.Control();
    
    if(not keyboardState['t']){
        camera.Move(dt);
        thunderbolt->MovePosition(dt);
        airscrew->MovePosition(dt);
    }

    thunderbolt->frenet(dt);

    glutPostRedisplay();
}

int main(int argc, char * argv[]) 
{
	glutInit(&argc, argv);
#if !defined(__APPLE__)
	glutInitContextVersion(majorVersion, minorVersion);
#endif
	glutInitWindowSize(windowWidth, windowHeight); 
	glutInitWindowPosition(50, 50);
#if defined(__APPLE__)
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_3_2_CORE_PROFILE);  
#else
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
#endif
	glutCreateWindow("3D Mesh Rendering");

#if !defined(__APPLE__)
	glewExperimental = true;	
	glewInit();
#endif
	printf("GL Vendor    : %s\n", glGetString(GL_VENDOR));
	printf("GL Renderer  : %s\n", glGetString(GL_RENDERER));
	printf("GL Version (string)  : %s\n", glGetString(GL_VERSION));
	glGetIntegerv(GL_MAJOR_VERSION, &majorVersion);
	glGetIntegerv(GL_MINOR_VERSION, &minorVersion);
	printf("GL Version (integer) : %d.%d\n", majorVersion, minorVersion);
	printf("GLSL Version : %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));
	
	onInitialization();

	glutDisplayFunc(onDisplay); 
	glutIdleFunc(onIdle);
	glutKeyboardFunc(onKeyboard);
	glutKeyboardUpFunc(onKeyboardUp);
	glutReshapeFunc(onReshape);

	glutMainLoop();
	onExit();
	return 1;
}

