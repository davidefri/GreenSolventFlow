// JS/scroll.js
window.onscroll = function() {scrollFunction()};

function scrollFunction() {
    var scrollToTopBtn = document.getElementById("scrollToTopBtn");
    if (document.body.scrollTop > 20 || document.documentElement.scrollTop > 20) {
        scrollToTopBtn.style.display = "block";
    } else {
        scrollToTopBtn.style.display = "none";
    }
}

// Quando l'utente clicca sul pulsante, scorri fino all'inizio del documento
document.getElementById("scrollToTopBtn").addEventListener("click", function() {
    document.body.scrollTop = 0; // Per Safari
    document.documentElement.scrollTop = 0; // Per Chrome, Firefox, IE e Opera
});