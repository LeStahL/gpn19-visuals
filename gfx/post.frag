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

void rand(in vec2 x, out float n)
{
    x += 400.;
    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec4 col = vec4(0.);
    float bound = sqrt(iFSAA)-1.;
    
    float delta;
    if(iEffect == 1) delta = .02;
    else delta = .0;
    
    float n;
    rand(floor(20.*fragCoord.y/iResolution.y*c.xx-1337.*floor(12.*iTime)),n);
    
    
   	for(float i = -.5*bound; i<=.5*bound; i+=1.)
        for(float j=-.5*bound; j<=.5*bound; j+=1.)
        {
            vec3 cl = texture(iChannel0, fragCoord/iResolution.xy+delta*n*c.xy+vec2(i,j)*3.0/max(bound,1.)/iResolution.xy).rgb,
                cr = texture(iChannel0, fragCoord/iResolution.xy-delta*n*c.xy+vec2(i,j)*3.0/max(bound,1.)/iResolution.xy).rgb,
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
