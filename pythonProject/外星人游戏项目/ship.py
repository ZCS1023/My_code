import pygame
from pygame.sprite import Sprite

class Ship(Sprite):

    def __init__(self, ai_settings, screen):

        #初始化飞船
        super(Ship, self).__init__()
        '''Python中的super(Name, self).__init__()是指首先找到Name的父类（比如是类NName），
        然后把类Name的对象self转换为类NName的对象，然后“被转换”的类NName对象调用自己的init函数，
        其实简单理解就是子类Name把父类的__init__()放到自己的__init__()当中，
        这样子类就有了父类的__init__()的各种元素'''
        self.screen = screen
        self.ai_settings = ai_settings

        #加载图片,并获取位置
        self.image = pygame.image.load('images/ship.bmp')
        self.rect = self.image.get_rect()
        self.screen_rect = screen.get_rect()

        self.rect.centerx = self.screen_rect.centerx
        self.rect.bottom = self.screen_rect.bottom

        self.center = float(self.rect.centerx)

        self.moving_right = False
        self.moving_left = False
        
    def center_ship(self):
        #使其居中
        self.center = self.screen_rect.centerx
        
    def update(self):

        if self.moving_right and self.rect.right < self.screen_rect.right:
            self.center += self.ai_settings.ship_speed_factor
        if self.moving_left and self.rect.left > 0:
            self.center -= self.ai_settings.ship_speed_factor

        self.rect.centerx = self.center

    def blitme(self):

        self.screen.blit(self.image, self.rect)
