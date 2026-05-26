import pygame.font
from pygame.sprite import Group
from ship import Ship

class Scoreboard():

    def __init__(self, ai_settings, screen, stats):

       #初始化相关属性

        self.screen = screen  #初始化屏幕
        self.screen_rect = screen.get_rect()  #获取屏幕的位置传递给screen_rect
        self.ai_settings = ai_settings   #初始化相关设置
        self.stats = stats   #初始化游戏统计
        
        #显示得分信息相关的字体设置
        self.text_color = (30, 30, 30)  #字体颜色
        self.font = pygame.font.SysFont('timesnewroman', 38,italic=True) #字体字型，大小等相关设置

        #准备得分图像(包含当前分数，最高得分，当前级别和当前飞船数)
        self.prep_score()
        self.prep_high_score()
        self.prep_level()
        self.prep_ships()

    def prep_score(self):
       #将得分(字符串)转换为图像(图)

        score_str ="Current score:"+str(int(round(self.stats.score, -1))) #将得分以离当前数字最近的10的倍数显示
        self.score_image = self.font.render(score_str, True, self.text_color,
            self.ai_settings.bg_color)
            
        #将得分放在屏幕的右上角
        self.score_rect = self.score_image.get_rect()
        self.score_rect.right = self.screen_rect.right - 20
        self.score_rect.top = 20
        
    def prep_high_score(self):
        #将最高得分转换为图像

        high_score_str = "Highest score:"+str(int(round(self.stats.high_score, -1)))
        self.high_score_image = self.font.render(high_score_str, True,
            self.text_color, self.ai_settings.bg_color)
                
        #将最高得分放置在屏幕顶部中间
        self.high_score_rect = self.high_score_image.get_rect()
        self.high_score_rect.centerx = self.screen_rect.centerx
        self.high_score_rect.top = self.score_rect.top
        
    def prep_level(self):
        #将当前级别转换为图像

        self.level_image = self.font.render(str(self.stats.level), True,
                self.text_color, self.ai_settings.bg_color)
        
        #将级别放在当前得分的下面
        self.level_rect = self.level_image.get_rect()
        self.level_rect.right = self.score_rect.right
        self.level_rect.top = self.score_rect.bottom + 10
        
    def prep_ships(self):
        #展示还剩几艘飞船
        self.ships = Group()
        for ship_number in range(self.stats.ships_left):
            ship = Ship(self.ai_settings, self.screen)  #使用Ship类创建ship实例
            ship.rect.x = 10 + ship_number * ship.rect.width
            ship.rect.y = 10   #设置ship位置
            self.ships.add(ship) #将ship实例加入到ships里
        
    def show_score(self):
        #将得分显示出来

        self.screen.blit(self.score_image, self.score_rect)  #显示当前得分
        self.screen.blit(self.high_score_image, self.high_score_rect)  #显示最高得分
        self.screen.blit(self.level_image, self.level_rect) #显示当前级别
        #在屏幕上画出飞船
        self.ships.draw(self.screen)
