$(document).ready(function() {

    var count = null;
    var counter = null;
    var highScore = 0;
    var currentScore = null;
    var gameActive = false;
    operator = document.getElementById("operator").value
    f = new Function('x', 'y', 'return x' + operator + 'y');

    document.getElementById("answer").disabled = true;

    function getParameterByName(name) {
        var match = RegExp('[?&]' + name + '=([^&]*)').exec(window.location.search);
        return match && decodeURIComponent(match[1].replace(/\+/g, ' '));
    }

    function timer() {
        count = count - 1;
        if (count < 0) {
            endGame();
            return;
        }

        document.getElementById("countdown").innerHTML = count;
    }

    function startGame() {

        $("#results").hide();
        $("#results").html("");
        gameActive = true;
        currentScore = 0;
        updateScore(currentScore);

        createQuestion();

        clearInterval(counter);

        count = 30;

        if (getParameterByName('time') != null) {
            count = getParameterByName('time');

        }

        document.getElementById("answer").disabled = false;
        document.getElementById("answer").value = "";
        document.getElementById("countdown").innerHTML = count;
        counter = setInterval(timer, 1000);

        document.getElementById("answer").focus();

    }

    function updateScore(newScore) {
        document.getElementById("currentScore").innerHTML =
            newScore;
    }

    var questionStartTime = null;

    function createQuestion() {

        do {

            //	        maximumOperand = 5;
            maximumOperand = 12;

            if (getParameterByName('operand') != null) {
                maximumOperand = getParameterByName('operand');
            }

            existing = document.getElementById("equation").innerHTML;

            // operand1 = Math.round(Math.random() * maximumOperand) * Math.pow(10, Math.round(Math.random()));
            // operand2 = Math.round(Math.random() * maximumOperand) * Math.pow(10, Math.round(Math.random()));
            operand1 = Math.round(Math.random() * maximumOperand);
            if (operator=='/' && operand1==0) {
                operand1=1;
            }
            operand2 = Math.round(Math.random() * maximumOperand);
            if (operator=='/' && operand2==0) {
                operand2=1;
            }


            if (getParameterByName('operator') != null) {
                operator = getParameterByName('operator');
                f = new Function('x', 'y', 'return x' + operator + 'y');
            }

            if (operator === '-' && operand2>operand1) {
                var placeholder = operand1;
                operand1 = operand2;
                operand2 = placeholder;
            }
            if (operator === '/') {
                var placeholder = operand1;
                operand1 = operand1*operand2;
                operand2 = placeholder;
            }

            if (getParameterByName('number') != null) {
                document.getElementById("equation").innerHTML = f(operand1, operand2).toString();
            } else {
                document.getElementById("equation").innerHTML = operand1 +
                    operator + operand2;
            }

            next = document.getElementById("equation").innerHTML;

            questionStartTime = new Date().getTime();

        } while (existing == next);

    }

    function evaluateAnswer() {
        if (gameActive) {

            var readyToEvaluate = f(operand1, operand2).toString().length == document.getElementById("answer").value.toString().length;

            if (readyToEvaluate) {

                if (parseInt(document.getElementById("answer").value) == f(operand1, operand2)) {

                    var correctSound = document.getElementById("correctSound");
                    correctSound.muted=false; // firefox warning
//                    correctSound.src = correctSound.src; // safari
                    correctSound.play();

                    currentScore++;
                    updateScore(currentScore);

                    var resultLog = "<div>";
                    resultLog += "<span class='infoHeader'>"+document.getElementById("equation").textContent+"&nbsp;=&nbsp;"+document.getElementById("answer").value+"</span>";
                    resultLog += "<span class='infoHeader'>Correct</span>";
                    resultLog += "</div>";

                    $("#results").append(resultLog);

                } else {

                    var incorrectSound = document.getElementById("incorrectSound");
                    incorrectSound.muted=false; // firefox warning
//                    incorrectSound.src = incorrectSound.src; // safari
                    incorrectSound.play();

                    var resultLog = "<div>";
                    resultLog += "<span class='infoHeader'>"+document.getElementById("equation").textContent+"&nbsp;!=&nbsp;"+document.getElementById("answer").value+"</span>";
                    resultLog += "<span class='infoHeader'>Incorrect</span>";
                    resultLog += "</div>";

                    $("#results").append(resultLog);

                }

                document.getElementById("answer").value = "";

                createQuestion();
                document.getElementById("answer").focus();

            }

        }

    }

    function endGame() {

        if (gameActive) {

            document.getElementById("countdown").innerHTML = 0;
            clearInterval(counter);
            gameActive = false;
            document.getElementById("answer").disabled = true;
            document.getElementById("equation").innerHTML = "Complete!";

            $("#results").show();

            if (currentScore > highScore) {

                var highScoreSound = document.getElementById("highScoreSound");
                highScoreSound.muted=false; // firefox warning
//                highScoreSound.src = highScoreSound.src; // safari
                highScoreSound.play();

                highScore = currentScore;
                document.getElementById("highScore").innerHTML =
                    highScore;

            } else {

                var gameOverSound = document.getElementById("gameOverSound");
                gameOverSound.muted=false; // firefox warning
//                gameOverSound.src = gameOverSound.src; // safari
                gameOverSound.play();

            }

            document.getElementById("answer").value = "";

        }

    }

    document.body.onkeyup = function(e) {
        if (e.keyCode == 32) {
            startGame();
        }
        if (e.keyCode == 27) {
            endGame();
        }
    };

    $("#answer").bind("keyup", function() {
        evaluateAnswer();
    });

    //DOM EVENTS
    $('#newGameBtn').on('click', function() {
        startGame();
    });

    $('#endGameBtn').on('click', function() {
        endGame();
    });
    $('#operator').on('change', function() {
        operator=document.getElementById("operator").value;
        f = new Function('x', 'y', 'return x' + operator + 'y');
        startGame();
    });

});
