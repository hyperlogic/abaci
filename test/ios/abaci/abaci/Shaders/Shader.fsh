//
//  Shader.fsh
//  abaci
//
//  Created by anthonythibault on 11/19/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

varying lowp vec4 colorVarying;

void main()
{
    gl_FragColor = colorVarying;
}
