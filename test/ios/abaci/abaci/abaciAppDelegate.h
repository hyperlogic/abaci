//
//  abaciAppDelegate.h
//  abaci
//
//  Created by anthonythibault on 11/19/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#import <UIKit/UIKit.h>

@class abaciViewController;

@interface abaciAppDelegate : NSObject <UIApplicationDelegate> {

}

@property (nonatomic, retain) IBOutlet UIWindow *window;

@property (nonatomic, retain) IBOutlet abaciViewController *viewController;

@end
