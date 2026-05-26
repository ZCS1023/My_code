import pygame
from pygame.sprite import Group

from settings import Settings
from game_stats import GameStats
from scoreboard import Scoreboard
from button import Button
from ship import Ship
import game_functions as gf

def run_game():

    #初始化游戏和设置，并创建一个屏幕对象
    pygame.init()
    ai_settings = Settings()
    screen = pygame.display.set_mode(
        (ai_settings.screen_width, ai_settings.screen_height))
    pygame.display.set_caption("Alien Invasion")
    
    #制作play按钮
    play_button = Button(ai_settings, screen, "Play")
    
    #创建实例，并创建记分牌
    stats = GameStats(ai_settings)
    sb = Scoreboard(ai_settings, screen, stats)
    
    #设置背景颜色
    bg_color = (230, 230, 230)
    
    #创建一艘船,一组子弹，和一组外星人
    ship = Ship(ai_settings, screen)
    bullets = Group()
    aliens = Group()
    
    #创建外星人群
    gf.create_fleet(ai_settings, screen, ship, aliens)

    #开始游戏的主循环
    while True:
        gf.check_events(ai_settings, screen, stats, sb, play_button, ship,
            aliens, bullets)  #不论何时，点击play按钮都应当有响应

        if stats.game_active:  #判断game_active指示器是否为True
            ship.update()
            gf.update_bullets(ai_settings, screen, stats, sb, ship, aliens,
                bullets)
            gf.update_aliens(ai_settings, screen, stats, sb, ship, aliens,
                bullets)  #更新外星人的同时将实参sb传递过去
        
        gf.update_screen(ai_settings, screen, stats, sb, ship, aliens,
            bullets, play_button)

run_game()
