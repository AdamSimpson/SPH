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

#ifndef CAMERA_H
#define CAMERA_H

#include "glm/glm.hpp"
#include "gl.hpp"

class Camera
{
    public:
        Camera(GL& gl); // Initializer list in .cpp
        // Binding index for global uniform matrices
        static const int g_GlobalMatricesBindingIndex = 0;
        static GLuint g_GlobalMatricesUBO;

        // Binding index for global uniform lights
        static const int g_GlobalLightBindingIndex = 1;
        static GLuint g_GlobalLightUBO;

        void rotate_yaw(float degrees);
        void rotate_pitch(float degrees);
        void rotate_yaw_pitch(float degrees_pitch, float degrees_yaw);
        void zoom_view(float d_zoom_factor);
        void move_along_view(float dx);
        void update_view();

    private:
        GL& gl;
        // Max degrees user allowed to rotate
        float max_degrees_rotate;
        // zoom factor for perspective matrix
        float zoom_factor;
        // Camera vectors
        glm::vec3 position;
        glm::vec3 position_default;
        glm::vec3 look_at;
};

#endif
