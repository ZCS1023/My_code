import math

def draw_moon(radius, message):
    """
    在控制台绘制一个圆月并显示祝福语
    """
    diameter = radius * 2

    # 计算中文字符宽度（每个汉字占2单位）
    def visual_length(s):
        return sum(2 if '\u4e00' <= c <= '\u9fff' else 1 for c in s)

    vis_len = visual_length(message)
    left_pad = (diameter - vis_len) // 2
    right_pad = diameter - vis_len - left_pad

    # 打印月亮
    for y in range(-radius, radius + 1):
        line = ""
        for x in range(-radius, radius + 1):
            adjusted_x = x * 0.9  # 微调宽高比
            distance = math.sqrt(adjusted_x**2 + y**2)
            if distance <= radius:
                if distance > radius - 1.2:
                    line += "*"
                else:
                    line += "●"  # 用实心点填充月亮内部
            else:
                line += " "
        print(line)

    # 打印祝福语（居中）
    print(" " * left_pad + message + " " * right_pad)

# 主程序
if __name__ == "__main__":
    draw_moon(12, "中秋快乐！")