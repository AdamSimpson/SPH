#ifndef GL_SHADER_FONTSTASH_H
#define GL_SHADER_FONTSTASH_H

#include "ogl_utils.h"

struct FONScontext* gl_shader_fonsCreate(int width, int height, int screen_width, int screen_height, int flags);
void gl_shader_fonsDelete(struct FONScontext* ctx);

unsigned int glfonsRGBA(unsigned char r, unsigned char g, unsigned char b, unsigned char a);

#endif

#ifdef GL_SHADER_FONTSTASH_IMPLEMENTATION

struct GLFONScontext {
    // Program handle
    GLuint program;

    // Locations
    GLint position_location;
    GLint tex_coord_location;
    GLint tex_location;
    GLint screen_dims_location;
    GLint color_location;

    // Required to have seperate buffers for verticies and tex coords
    GLuint vbo_vert, vbo_tex, vbo_color;

    // Uniforms
    GLuint tex;
    GLuint screen_dims;

    // ATLAS width and height
    int width, height;

    // Screen pixel dimensions
    int screen_width, screen_height;
};

static void glfons__create_shaders(void * userPtr)
{
    struct GLFONScontext* gl = (struct GLFONScontext*)userPtr;

    // Compile vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    #ifdef GLES
      compile_shader(vertexShader, "SPH/font_es.vert");
    #else
      compile_shader(vertexShader, "font.vert");
    #endif

    // Compile frag shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    #ifdef GLES
      compile_shader(fragmentShader, "SPH/font_es.frag");
    #else
      compile_shader(fragmentShader, "font.frag");
    #endif

    // Create shader program
    gl->program = glCreateProgram();
    glAttachShader(gl->program, vertexShader);
    glAttachShader(gl->program, fragmentShader);
   
    // Link and use program
    glLinkProgram(gl->program);

    // Get position location
    gl->position_location = glGetAttribLocation(gl->program, "position");
    // Get tex_coord location
    gl->tex_coord_location = glGetAttribLocation(gl->program, "tex_coord");
    // Get color location
    gl->color_location = glGetAttribLocation(gl->program, "color");
    // Get tex uniform location
    gl->tex_location = glGetUniformLocation(gl->program, "tex");
    // Get screen dims uniform location
    gl->screen_dims_location = glGetUniformLocation(gl->program, "screen_dims");
}

static void glfons__create_buffers(void * userPtr)
{
    struct GLFONScontext* gl = (struct GLFONScontext*)userPtr;

    // Must use VAO with VBO
    // VAO is REQUIRED for OpenGL 3+ when using VBO I believe
    #ifndef GLES
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    #endif

    glGenBuffers(1, &gl->vbo_vert);
    glGenBuffers(1, &gl->vbo_tex);
    glGenBuffers(1, &gl->vbo_color);
}

static int glfons__renderCreate(void* userPtr, int width, int height)
{
	struct GLFONScontext* gl = (struct GLFONScontext*)userPtr;
	glGenTextures(1, &gl->tex);
	if (!gl->tex) return 0;
	gl->width = width;
	gl->height = width;
	glBindTexture(GL_TEXTURE_2D, gl->tex);
    #ifdef GLES
 	glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, gl->width, gl->height, 0, GL_ALPHA, GL_UNSIGNED_BYTE, 0);
	#else
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, gl->width, gl->height, 0, GL_RED, GL_UNSIGNED_BYTE, 0);
    #endif

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	// Programable pipeline requires explicit buffer and shader managmement
        glfons__create_buffers(gl);
        glfons__create_shaders(gl);

	return 1;
}

static void glfons__renderUpdate(void* userPtr, int* rect, const unsigned char* data)
{
	struct GLFONScontext* gl = (struct GLFONScontext*)userPtr;
	int w = rect[2] - rect[0];
	int h = rect[3] - rect[1];

	if (gl->tex == 0) return;
	glBindTexture(GL_TEXTURE_2D, gl->tex);
	glPixelStorei(GL_UNPACK_ALIGNMENT,1);

	#ifdef GLES
    int y;
    for(y=0; y<h; y++) {
	    const unsigned char *row = data + (rect[1]+y)*gl->width + rect[0];
	    glTexSubImage2D(GL_TEXTURE_2D, 0, rect[0], rect[1]+y, w, 1, GL_ALPHA, GL_UNSIGNED_BYTE, row);
    }
	#else
	glPixelStorei(GL_UNPACK_ROW_LENGTH, gl->width);
	glPixelStorei(GL_UNPACK_SKIP_PIXELS, rect[0]);
	glPixelStorei(GL_UNPACK_SKIP_ROWS, rect[1]);
	glTexSubImage2D(GL_TEXTURE_2D, 0, rect[0], rect[1], w, h, GL_RED,GL_UNSIGNED_BYTE, data);
	#endif

}

static void glfons__renderDraw(void* userPtr, const float* verts, const float* tcoords, const unsigned int* colors, int nverts)
{
	struct GLFONScontext* gl = (struct GLFONScontext*)userPtr;
	if (gl->tex == 0) return;

        // Size of each vertex in bytes
        size_t vert_size = 2*sizeof(GL_FLOAT);

        // Draw texture
        glUseProgram(gl->program);

        // Bind vert buffer
        glBindBuffer(GL_ARRAY_BUFFER, gl->vbo_vert);
        // Fill vert buffer
        glBufferData(GL_ARRAY_BUFFER, nverts*vert_size, verts, GL_DYNAMIC_DRAW);
	    // Get and enable vertex pointer
        glVertexAttribPointer(gl->position_location, 2, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(gl->position_location);

        // Bind tex coord buffer
        glBindBuffer(GL_ARRAY_BUFFER, gl->vbo_tex);
        // Fill tex coord buffer
        glBufferData(GL_ARRAY_BUFFER, nverts*vert_size, tcoords, GL_DYNAMIC_DRAW);
  	    // Get and enable tex coord pointer
        glVertexAttribPointer(gl->tex_coord_location, 2, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(gl->tex_coord_location);

        // Bind color buffer
        glBindBuffer(GL_ARRAY_BUFFER, gl->vbo_color);
        // Fill color buffer
        glBufferData(GL_ARRAY_BUFFER, 4*nverts, colors, GL_DYNAMIC_DRAW);
        // Get and enable color pointer
        glVertexAttribPointer(gl->color_location, 4, GL_UNSIGNED_BYTE, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(gl->color_location);

 	    // Bind tex
        glBindTexture(GL_TEXTURE_2D, gl->tex);
	    //Set tex uniform
        glUniform1i(gl->tex_location, 0);

	    // Set screen dimension uniform
        glUniform2i(gl->screen_dims_location, gl->screen_width, gl->screen_height);

        glDrawArrays(GL_TRIANGLES, 0, nverts);

}

static void glfons__renderDelete(void* userPtr)
{
	struct GLFONScontext* gl = (struct GLFONScontext*)userPtr;
	if (gl->tex)
		glDeleteTextures(1, &gl->tex);
	gl->tex = 0;
	free(gl);

        // Need to cleanup here....
}

struct FONScontext* gl_shader_fonsCreate(int width, int height, int screen_width, int screen_height, int flags)
{
	struct FONSparams params;
	struct GLFONScontext* gl;

	gl = (struct GLFONScontext*)malloc(sizeof(struct GLFONScontext));
	if (gl == NULL) goto error;
	memset(gl, 0, sizeof(struct GLFONScontext));

	memset(&params, 0, sizeof(params));
	params.width = width;
	params.height = height;
	params.flags = flags;
	params.renderCreate = glfons__renderCreate;
	params.renderUpdate = glfons__renderUpdate;
	params.renderDraw = glfons__renderDraw; 
	params.renderDelete = glfons__renderDelete;
	params.userPtr = gl;

        gl->screen_width = screen_width;
        gl->screen_height = screen_height;

	return fonsCreateInternal(&params);

error:
	if (gl != NULL) free(gl);
	return NULL;
}

void gl_shader_fonsDelete(struct FONScontext* ctx)
{
	fonsDeleteInternal(ctx);
}

unsigned int glfonsRGBA(unsigned char r, unsigned char g, unsigned char b, unsigned char a)
{
	return (r) | (g << 8) | (b << 16) | (a << 24);
}

#endif
