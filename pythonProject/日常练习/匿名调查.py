class AnonySurvey():
    def __init__(self,question):
        self.question=question
        self.responses=[]
    def show_question(self):
        print(self.question)
    def store_response(self,new_response):
        self.responses.append(new_response)
    def show_response(self):
        print("\n收集到的结果为:")
        for response in self.responses:
            print(response)

'''my_survey=AnonySurvey("你的母语是？")
my_survey.show_question()
print("输入q来退出程序")
while True:
    response=input("母语：")
    if response=='q':
        break
    my_survey.store_response(response)

my_survey.show_response()'''