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

uniform float iFader0;
uniform float iFader1;
uniform float iFader2;
uniform float iFader3;
uniform float iFader4;
uniform float iFader5;
uniform float iFader6;
uniform float iFader7;

uniform float iDial0;
uniform float iDial1;
uniform float iDial2;
uniform float iDial3;
uniform float iDial4;
uniform float iDial5;
uniform float iDial6;
uniform float iDial7;


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

void dvoronoi(in vec2 x, out float d, out vec2 z)
{
    vec2 y = floor(x);
       float ret = 1.;
    vec2 pf=c.yy, p;
    float df=10.;
    
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            float pa;
            rand(p, pa);
            p += pa;
            
            d = length(x-p);
            
            if(d < df)
            {
                df = d;
                pf = p;
            }
        }
    for(int i=-1; i<=1; i+=1)
        for(int j=-1; j<=1; j+=1)
        {
            p = y + vec2(float(i), float(j));
            float pa;
            rand(p, pa);
            p += pa;
            
            vec2 o = p - pf;
            d = length(.5*o-dot(x-pf, o)/dot(o,o)*o);
            ret = min(ret, d);
        }
    
    d = ret;
    z = pf;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec4 col = vec4(0.);
    float bound = sqrt(iFSAA)-1.;
    
    float delta = 0.;
    vec2 n;
    
    // Chromatic distortion
    if(iFader0 > 0.) 
    {
        delta = mix(.0,.02,iFader0);
        rand(floor(20.*fragCoord.y/iResolution.y*c.xx-1337.*floor(12.*iTime)),n.x);
        rand(floor(20.*fragCoord.y/iResolution.y*c.xx-1337.*floor(12.*iTime)+2337.),n.y);
    }
    
    // HF noise
    if(iFader1 > 0.)
    {
        lfnoise(12.*fragCoord-iTime, n.x);
        lfnoise(12.*fragCoord-iTime-1337., n.y);
        fragCoord += mix(1.,20.,iFader1)*n;
    }
    
    // LF noise
    if(iFader2 > 0.)
    {
        lfnoise(22.*fragCoord/iResolution-3.*iTime, n.x);
        lfnoise(22.*fragCoord/iResolution-3.*iTime-1337., n.y);
        fragCoord += mix(0.,22.,iFader2)*n;
    }
    
    // Kaleidoscope
    if(iFader3 > 0.)
    {
        float a = iResolution.x/iResolution.y;
        vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
//         rand(floor(.33*iTime)*c.xx, n.x);
//         n.x = max(floor(12.*n.x),3.);
        n.x = floor(mix(3.,10.,iFader3));
        float phi = abs(mod(atan(uv.y, uv.x),pi/n.x)-.5*pi/n.x);
        uv = length(uv)*vec2(cos(phi), sin(phi));
        fragCoord = (uv + .5*vec2(a,1.))*iResolution.yy;
    }
    
    if(iFader4 > 0.)
    {
        float a = iResolution.x/iResolution.y;
        vec2 uv = fragCoord/iResolution.yy-0.5*vec2(a, 1.0);
        
        float dv;
        vec2 ind;
        dvoronoi(mix(1.,100.,1.-iFader4)*uv, dv, ind);
        uv = ind/mix(1.,100.,1.-iFader4);
        
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
