/*
The MIT License (MIT)

Copyright (c) 2014 Adam Simpson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#ifndef GLFW_UTILS_H
#define GLFW_UTILS_H

extern "C"
{
    #include <GL/glew.h>
    #include <GLFW/glfw3.h>
}

#include "structs.h"

class GLFW
{
    public:
        GLFW(): cursor_x(0), cursor_y(0), cursor_view_x(0), cursor_view_y(0);
        ~GLFW();
        void check_user_input();
        void swap_buffers();
        bool window_should_close();
        void exit_program();
        void pixel_to_gl(const int pixel_x, const int pixel_y, float &gl_x, float &gl_y);
        int screen_width() const { return this->screen_width };
        int screen_height() const { return this->screen_height};
    private:
        GLFWwindow *window;
        float screen_width;
        float screen_height;
        // GLFW mouse coordinates for "normal" cursor coordinates
        float cursor_x;
        float cursor_y;
        // GLFW mouse coordinates for "view" mode cursor coordinates
        float cursor_view_x;
        float cursor_view_y;
}

// GLFW callback function prototypes
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
static void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void wheel_callback(GLFWwindow* window, double x, double y);
void error_callback(int error, const char* description);

#endif
