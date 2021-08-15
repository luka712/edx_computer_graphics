# version 330 core
// Do not use any version older than 330!

/* This is the fragment shader for reading in a scene description, including 
   lighting.  Uniform lights are specified from the main program, and used in 
   the shader.  As well as the material parameters of the object.  */

// Inputs to the fragment shader are the outputs of the same name of the vertex shader.
// Note that the default output, gl_Position, is inaccessible!
in vec3 mynormal; 
in vec4 myvertex; 

// You will certainly need this matrix for your lighting calculations
uniform mat4 modelview;

// This first defined output of type vec4 will be the fragment color
out vec4 fragColor;

uniform vec3 color;

const int numLights = 10; 
uniform bool enablelighting; // are we lighting at all (global).
uniform vec4 lightposn[numLights]; // positions of lights 
uniform vec4 lightcolor[numLights]; // colors of lights
uniform int numused;               // number of lights used

// Now, set the material parameters.
// I use ambient, diffuse, specular, shininess. 
// But, the ambient is just additive and doesn't multiply the lights.  

uniform vec4 ambient; 
uniform vec4 diffuse; 
uniform vec4 specular; 
uniform vec4 emission; 
uniform float shininess; 

void main (void) 
{       
    if (enablelighting) {       
        vec4 finalcolor; 

        // YOUR CODE FOR HW 2 HERE
        // A key part is implementation of the fragment shader

        vec4 diffuse_sum = vec4(0.0);
        vec4 point_sum = vec4(0.0);
       
        // TODO: note that some code here can be optimized, ie normal should be normalized earlier.
        vec4 vertex = modelview * myvertex;
        vec4 normal = transpose(inverse(modelview)) * vec4(mynormal, 1.0);
        for(int i = 0; i < numLights;i++)
        { 
            // get the light position and init light dir.
            vec3 light_pos = (lightposn[i]).xyz;
            vec3 light_dir = vec3(0.0);
            
            // if w > 0.0 use calc for point light
            if(lightposn[i].w > 0.0) 
            {
                // get vector pointing to light
                light_dir = normalize(light_pos - vertex.xyz);
                // amount of light between normal and light dir.
                float n_dot_l = dot(normalize(normal.xyz), light_dir);
                // must be at least 0.0, no negative light
                n_dot_l = max(n_dot_l, 0.0);
                diffuse_sum += diffuse * lightcolor[i] * n_dot_l;
            }
            else // directional light
            {
                // light dir is just light position in case of directional light.
                light_dir = normalize(light_pos.xyz);
                // amounf of light between normal and light dir
                float n_dot_l = dot(normalize(normal.xyz), light_dir);
                n_dot_l = max(n_dot_l, 0.0);
                diffuse_sum += diffuse * lightcolor[i] * n_dot_l;
            }
            
            
            // specular component goes here, it uses blinn phong with half vector between eye direction and light dir
            // https://en.wikipedia.org/wiki/Blinn%E2%80%93Phong_reflection_model
            const vec3 eye = vec3(0.0);
            // homogenize
            vec3 position = vertex.xyz / vertex.w;
            vec3 eye_dir = normalize(eye - position);
            vec3 half_v_normalized = normalize(eye_dir.xyz + light_dir);
            
            float n_dot_h = dot( normalize(normal.xyz), half_v_normalized);
            n_dot_h = max(n_dot_h, 0.0);
            point_sum += specular * lightcolor[i] * pow(n_dot_h, shininess);
        }

        // Color all pixels black for now, remove this in your implementation!
        finalcolor = ambient + emission + diffuse_sum + point_sum;

        fragColor = finalcolor; 
    } else {
        fragColor = vec4(color, 1.0f); 
    }
}
