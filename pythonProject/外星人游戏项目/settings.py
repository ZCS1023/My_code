class Settings():
    #存储alien_invasion的所有设置

    def __init__(self):

        #屏幕设置
        self.screen_width = 1200
        self.screen_height = 800
        self.bg_color = (230, 230, 230)
        
        #飞船设置
        self.ship_limit = 3
            
        #子弹设置
        self.bullet_width = 3
        self.bullet_height = 15
        self.bullet_color = 255, 0, 0
        self.bullets_allowed = 3
        
        #外星人下落速度设置
        self.fleet_drop_speed = 10
            
        #游戏速度提升多少
        self.speedup_scale = 1.1
        #外星人射杀点数提升多少
        self.score_scale = 1.5
    
        self.initialize_dynamic_settings()

    def initialize_dynamic_settings(self):
        #初始化游戏全程变化的设置

        self.ship_speed_factor = 1.5
        self.bullet_speed_factor = 3
        self.alien_speed_factor = 0.2
        
        #得分
        self.alien_points = 50
    
        #1代表向右，-1向左
        self.fleet_direction = 1
        
    def increase_speed(self):
        #提高速度设置和外星人的点数
        self.ship_speed_factor *= self.speedup_scale
        self.bullet_speed_factor *= self.speedup_scale
        self.alien_speed_factor *= self.speedup_scale
        
        self.alien_points = int(self.alien_points * self.score_scale)
