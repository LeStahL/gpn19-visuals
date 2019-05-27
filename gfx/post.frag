/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19
 * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
 
#version 130

uniform float iFSAA;
uniform vec2 iResolution;
uniform float iTime;
uniform sampler2D iChannel0;
uniform int iEffect;

out vec4 gl_FragColor;

const vec3 c = vec3(1.,0.,-1.);
const float pi = acos(-1.);

void rand(in vec2 x, out float n)
{
    x += 400.;
    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);
}

void lfnoise(in vec2 t, out float n)
{
    vec2 i = floor(t);
    t = fract(t);
    t = smoothstep(c.yy, c.xx, t);
    vec2 v1, v2;
    rand(i, v1.x);
    rand(i+c.xy, v1.y);
    rand(i+c.yx, v2.x);
    rand(i+c.xx, v2.y);
    v1 = c.zz+2.*mix(v1, v2, t.y);
    n = mix(v1.x, v1.y, t.x);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec4 col = vec4(0.);
    float bound = sqrt(iFSAA)-1.;
    
    float delta;
    vec2 n;
    
    // Chromatic distortion
    if(iEffect == 1) 
    {
        delta = .02;
        rand(floor(20.*fragCoord.y/iResolution.y*c.xx-1337.*floor(12.*iTime)),n.x);
        rand(floor(20.*fragCoord.y/iResolution.y*c.xx-1337.*floor(12.*iTime)+2337.),n.y);
    }
    else delta = .0;
    
    // HF noise
    if(iEffect == 2)
    {
        lfnoise(12.*fragCoord-iTime, n.x);
        lfnoise(12.*fragCoord-iTime-1337., n.y);
        fragCoord += 20.*n;
    }
    
    // LF noise
    else if(iEffect == 3)
    {
        lfnoise(22.*fragCoord/iResolution-3.*iTime, n.x);
        lfnoise(22.*fragCoord/iResolution-3.*iTime-1337., n.y);
        fragCoord += 22.*n;
    }
    
    // Kaleidoscope
    else if(iEffect == 4)
    {
        float a = iResolution.x/iResolution.y;
        vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
        rand(floor(.33*iTime)*c.xx, n.x);
        n.x = max(floor(12.*n.x),3.);
        float phi = abs(mod(atan(uv.y, uv.x),pi/n.x)-.5*pi/n.x);
        uv = length(uv)*vec2(cos(phi), sin(phi));
        fragCoord = (uv + .5*vec2(a,1.))*iResolution.yy;
    }
    
   	for(float i = -.5*bound; i<=.5*bound; i+=1.)
        for(float j=-.5*bound; j<=.5*bound; j+=1.)
        {
            vec3 cl = texture(iChannel0, fragCoord/iResolution.xy+delta*n+vec2(i,j)*3.0/max(bound,1.)/iResolution.xy).rgb,
                cr = texture(iChannel0, fragCoord/iResolution.xy-delta*n+vec2(i,j)*3.0/max(bound,1.)/iResolution.xy).rgb,
                cc = texture(iChannel0, fragCoord/iResolution.xy+vec2(i,j)*3.0/max(bound,1.)/iResolution.xy).rgb;
            col += vec4(cl.r, cc.g, cr.b,1.);
        }
    col /= iFSAA;
    fragColor = col;
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
