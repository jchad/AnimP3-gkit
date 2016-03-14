
#version 330

#ifdef VERTEX_SHADER
in vec3 position;
in vec2 texcoord;

uniform mat4 mvpMatrix;

out vec2 vertex_texcoord;

void main( )
{
    gl_Position= mvpMatrix * vec4(position, 1);
    vertex_texcoord= texcoord;
}
#endif


#ifdef FRAGMENT_SHADER

in vec2 vertex_texcoord;

uniform sampler2D texture0;

void main( )
{
    vec4 color= texture(texture0, vertex_texcoord, 0);
    
    gl_FragColor= color;
}
#endif
