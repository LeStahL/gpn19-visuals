/* Generated with shader-compressor by NR4/Team210. */
#ifndef POST_H
#define POST_H
const char * post_frag =
"/* Endeavor by Team210 - 64k intro by Team210 at Revision 2k19\n"
" * Copyright (C) 2019  Alexander Kraus <nr4@z10.info>\n"
" *\n"
" * This program is free software: you can redistribute it and/or modify\n"
" * it under the terms of the GNU General Public License as published by\n"
" * the Free Software Foundation, either version 3 of the License, or\n"
" * (at your option) any later version.\n"
" *\n"
" * This program is distributed in the hope that it will be useful,\n"
" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
" * GNU General Public License for more details.\n"
" *\n"
" * You should have received a copy of the GNU General Public License\n"
" * along with this program.  If not, see <https://www.gnu.org/licenses/>.\n"
" */\n"
" \n"
"#version 130\n"
"\n"
"uniform float iFSAA;\n"
"uniform vec2 iResolution;\n"
"uniform float iTime;\n"
"uniform sampler2D iChannel0;\n"
"uniform int iEffect;\n"
"\n"
"out vec4 gl_FragColor;\n"
"\n"
"const vec3 c = vec3(1.,0.,-1.);\n"
"\n"
"void rand(in vec2 x, out float n)\n"
"{\n"
"    x += 400.;\n"
"    n = fract(sin(dot(sign(x)*abs(x) ,vec2(12.9898,78.233)))*43758.5453);\n"
"}\n"
"\n"
"void mainImage( out vec4 fragColor, in vec2 fragCoord )\n"
"{\n"
"    vec4 col = vec4(0.);\n"
"    float bound = sqrt(iFSAA)-1.;\n"
"    \n"
"    float delta;\n"
"    if(iEffect == 1) delta = .02;\n"
"    else delta = .0;\n"
"    \n"
"    float n;\n"
"    rand(floor(20.*fragCoord.y/iResolution.y*c.xx-1337.*floor(12.*iTime)),n);\n"
"    \n"
"    \n"
"   	for(float i = -.5*bound; i<=.5*bound; i+=1.)\n"
"        for(float j=-.5*bound; j<=.5*bound; j+=1.)\n"
"        {\n"
"            vec3 cl = texture(iChannel0, fragCoord/iResolution.xy+delta*n*c.xy+vec2(i,j)*3.0/max(bound,1.)/iResolution.xy).rgb,\n"
"                cr = texture(iChannel0, fragCoord/iResolution.xy-delta*n*c.xy+vec2(i,j)*3.0/max(bound,1.)/iResolution.xy).rgb,\n"
"                cc = texture(iChannel0, fragCoord/iResolution.xy+vec2(i,j)*3.0/max(bound,1.)/iResolution.xy).rgb;\n"
"            col += vec4(cl.r, cc.g, cr.b,1.);\n"
"        }\n"
"    col /= iFSAA;\n"
"    fragColor = col;\n"
"}\n"
"\n"
"void main()\n"
"{\n"
"    mainImage(gl_FragColor, gl_FragCoord.xy);\n"
"}\n"
"\n"
;
#endif
