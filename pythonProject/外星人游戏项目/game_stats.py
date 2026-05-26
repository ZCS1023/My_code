class GameStats():
    #记录alien_invation中的统计量
    
    def __init__(self, ai_settings):
        #初始化统计量
        self.ai_settings = ai_settings
        self.reset_stats()
        
        #在game-active指示器为false时开始游戏
        self.game_active = False
        
        #最高分在任何情况下都不会被重置
        self.high_score = 0  #除了关闭程序，否则最高分并不会重置
        
    def reset_stats(self):
        #初始化在游戏中可以改变的统计量
        self.ships_left = self.ai_settings.ship_limit
        self.score = 0
        self.level = 1
