attribute vec2 position;

uniform float radius;
uniform vec2 center;
uniform vec3 color;

varying vec2 local_frag_coord;

void main() {

    switch(gl_VertexID)
    {
        case 0:
            //Bottom-left
            local_frag_coord = vec2(-1.0, -1.0);
            break;
        case 1:
            //Top-left
            local_frag_coord = vec2(-1.0, 1.0);
            break;
        case 2:
            //Bottom-right
            local_frag_coord = vec2(1.0, -1.0);
            break;
        case 3:
            //Top-right
            local_frag_coord = vec2(1.0, 1.0);
            break;
    }

    gl_Position = vec4(position, 0.0, 1.0);
}

